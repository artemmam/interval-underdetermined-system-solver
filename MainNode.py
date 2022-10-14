from mpi4py import MPI
import time
import numpy as np
from check_box import diam, separate_box

class _task_wrapper:
	""" A wrapper for result object """

	def __init__(self, function, item, args):
		self.function = function
		self.args = [item] + args
		self.item = item


class MainNode:
	""" MainNode module

	Main nodes
	"""

	def __init__(self, function, queue=[], args=[], function_init='',
				 function_analysis='', debug=False, eps_bnb = 1):
		""" Main node cycle

		Parameters
		==========
		function :
		   The main calculation function

		function_init :
		   Initialization function

		function_analysis :
		   A function executed when a result is executed.

		Variables (Results)
		===================
		queues :
		   A list of queues

		results:
		   The stored results

		errors:
		   The stored errors
		"""

		# task and queues
		self.function = function
		self.function_init = function_init
		self.queues = queue
		self.args = args

		# connection (MPI)
		self.comm = MPI.COMM_WORLD
		self.size = self.comm.Get_size()
		self.rank = self.comm.Get_rank()
		self.status = MPI.Status()

		# internal parameters
		self.control_signal = ''
		self.wait_time = 5
		self.debug = debug
		self.working_node = {}

		# result variable
		self.results = []
		self.errors = []
		self.inside_boxes = []
		self.border_boxes = []
		self.eps_bnb = eps_bnb
		# initialize
		if self.function_init != '':
			self.function_init(self)

	def execute(self):
		""" Main node cycle

		"""

		# main-node cycle
		while True:
			# receive sub-node status
			signal, src = self._receive_status_signal()
			if signal == 'ready':
				self._send_task(src)
			elif signal == 'done':
				self._receive_result()
			elif signal == 'error':
				self._error_treatment()
			if len(self.working_node) == 0 and len(self.queues) == 0:
				self.results = [self.inside_boxes] + [self.border_boxes]
				break

	def terminate(self):
		""" Termination code

		Check whether queues are empty and all jobs are finished
		"""
		while True:
			if self._finish_condition():
				break

	def _receive_status_signal(self):
		signal = self.comm.recv(source=MPI.ANY_SOURCE, tag=0,
								status=self.status)
		src = self.status.source
		if self.debug:
			print('[main] receive status signal %s from %s' % (signal, src))
		return (signal, src)

	# def _analyze(self, result):
	#     time.sleep(3)
	#     if result % 2 != 0:
	#         print(self.queues)
	#         buf = self.queues
	#         buf.append([np.random.randint(1, 5), np.random.randint(6, 10)])
	#         self.queues = buf
	#         print('[main] append to queue :', self.queues)
	#     else:
	#         print(self.queues)
	#         buf = self.queues
	#         buf.append([np.random.randint(11, 16), np.random.randint(17, 25)])
	#         self.queues = buf
	#         print('[main] append to queue :', self.queues)

	def _send_control_signal(self, dest, signal):
		self.comm.send(signal, dest=dest, tag=0)

	def _send_task(self, src):
		if self.queues == []:
			print('[main] empty queue')
			self._send_control_signal(src, 'wait')
		else:
			item = self.queues[0]
			self.queues = self.queues[1:]
			args = self.args
			self._send_control_signal(src, 'task')
			self.comm.send(_task_wrapper(self.function, item, args), dest=src, tag=1)
			self.working_node[src] = time.time()

	def _receive_result(self):
		result = self.comm.recv(source=self.status.Get_source(), tag=2,
								status=self.status)
		src = self.status.source
		self.working_node.pop(src)
		if self.debug:
			print('[main] Result from node %s: %s' % (src, result.result))
		if result.result == "inside":
			self.inside_boxes.append(result.item)
		elif result.result == "border":
			if diam(result.item)<self.eps_bnb:
				self.border_boxes.append(result.item)
			else:
				box_l, box_r = separate_box(result.item)
				self.queues.append(box_l)
				self.queues.append(box_r)

	# if result.result == "inside":
	# self.inside_boxes.append(result.item)
	# elif result.result == "border":
	# if diam(result.item) < self.eps_bnb:
	# self.border_boxes.append(result.item)
	# else:

	# self._analyze(result.result)
	# else:
	# self.results.append(result.result)

	def _error_treatment(self):
		task = self.comm.recv(source=self.status.Get_source(), tag=2,
							  status=self.status)
		src = self.status.source
		self.working_node.pop(src)
		self.errors.append(task.args)

	def _finish_condition(self):
		if self.queues != []:
			return False
		if len(self.working_node) != 0:
			return False

		if self.debug:
			print('[main] All tasks finished, execute termination code')

		for i in range(self.size):
			if i != 0:
				if self.debug:
					print('[main] Send termination signal to node %d' % i)
				self._send_control_signal(i, 'end')
		return True