from scipy.stats import pointbiserialr
from argparse import ArgumentParser
from analyse_results import parse_all_results
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
		x = [r.tampered for r in results]
		y = [self.get_val(r) for r in results]

		print(x, y)
		brk()

		pb = pointbiserialr(x, y)
		return pb.correlation, pb.pvalue


class MutsPer(Correlation):
	def get_val(self, r):
		if r.total_sites:
			return float(r.muts_in_sites) / r.total_sites
		else:
			return 0


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)

	args = ap.parse_args()
	results = parse_all_results(args.fname[0])

	correlations = [
			Correlation("muts_in_sites"),
			Correlation("total_sites"),
			Correlation("total_singles"),
			MutsPer("muts_per")
			]


	for k, v in results.items():
		print(k)
		for c in correlations:
			print(c.name, c.calc(v))


if __name__ == "__main__":
	main()
