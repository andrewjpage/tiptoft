import re
class Gene:
	def __init__(self,name, kmers_with_coverage, kmers_without_coverage):
		self.name = name
		self.kmers_with_coverage = kmers_with_coverage
		self.kmers_without_coverage = kmers_without_coverage
		
	def __str__(self):
		s = self.name
		if not self.is_full_coverage():
			s += '*'
		s += "\t"+str(self.percentage_coverage())
		return s
		
	def is_full_coverage(self):
		if self.kmers_without_coverage  == 0:
			return True
		else:
			return False
			
	def percentage_coverage(self):
		total_kmers = self.kmers_with_coverage + self.kmers_without_coverage
		if total_kmers > 0:
			return int((self.kmers_with_coverage*100)/total_kmers)
		else:
			return 0
			