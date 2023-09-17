from scipy.stats import pointbiserialr, binom_test
from argparse import ArgumentParser
from analyse_results import parse_all_results
from math import sqrt
from pdb import set_trace as brk


class Correlation:
	def __init__(self, name, fn=None):
		self.name = name
		self.fn = fn

	def get_val(self, result):
		if self.fn:
			return self.fn(result)
		else:
			return getattr(result, self.name)

	def calc(self, results):
		x, y = [], []

		for r in results:
			val = self.get_val(r)
			if val is not None:
				x.append(r.tampered)
				y.append(val)

		pb = pointbiserialr(x, y)
		return pb.correlation, pb.pvalue

	def classify(self, results):
		values = [self.get_val(r) for r in results]
		values = [v for v in values if v is not None]

		# How well would this work as a classifier if we looked at whether you
		# were above or below the mean?
		total = sum(values)
		mean = float(total) / len(values)

		false_pos, true_pos = 0, 0
		false_neg, true_neg = 0, 0
		total_pos, total_neg = 0, 0

		# Now figure out the sensitvity and specificity
		for r in results:
			val = self.get_val(r)
			if val is None: continue
			estimation = val > mean
			truth = r.tampered

			if truth:
				total_pos += 1
			else:
				total_neg += 1

			if estimation == truth:
				if truth:
					true_pos += 1
				else:
					true_neg += 1
			else:
				if truth:
					false_neg += 1
				else:
					false_pos += 1

		sens = float(true_pos) / total_pos
		spec = float(true_neg) / total_neg
		return mean, sens, spec


class MutsPer(Correlation):
	def get_val(self, r):
		if r.total_sites:
			return float(r.muts_in_sites) / r.total_sites
		else:
			return None

def split_references(results):
	references, new_results = {}, {}
	for k, v in results.items():
		if k.startswith("WH1-"):
			assert len(v) == 1
			references[k] = v[0]
		else:
			new_results[k] = v

	return references, new_results


def rank(references, results, field):
	"""Considering genome, where do the results we found rank among the false
	simulated results?"""
	total, total_more = 0, 0
	for k, result in references.items():
		ref_val = getattr(result, field)
		simulated = results[k[4:]]	# Just remove the WH-1 at the start

		for result in simulated:
			if not result.tampered:
				continue
			total += 1
			if getattr(result, field) > ref_val:
				total_more += 1

		yield k, total, total_more


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)

	args = ap.parse_args()
	results = parse_all_results(args.fname[0])
	references, results = split_references(results)

	correlations = [
			Correlation("muts_in_sites"),
			Correlation("total_sites"),
			Correlation("total_singles"),
# 			MutsPer("muts_per")
			]


	for k, v in results.items():
		print(k)
		for c in correlations:
			print("{}: {:.4f} {:.4g}".format(c.name, *c.calc(v)))

			mean, sens, spec = c.classify(v)
			print("Mean: {:2f} Sensitivity: {:.2f} Specificity {:.2f}".format(
				mean, sens*100, spec*100))

			for k, ref in references.items():
				val = getattr(ref, c.name)
				print("{}: {} {}".format(k, val, val > mean))

	print("""
Where the real alignments for non-tampered genomes rank on each scores compared
to the simulated ones. We're counting how many non-tampered alignments score
higher:
""")
	for field in ("muts_in_sites", "total_sites", "total_singles"):
	# for field in ("total_sites",):
		for name, total, total_more in rank(references, results, field):
			print("{} {}: {}/{} {:.4f}". format(name, field, total_more,
				total, float(total_more) / total))


if __name__ == "__main__":
	main()
