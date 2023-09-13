from scipy.stats import pointbiserialr
from collections import namedtuple, defaultdict
from argparse import ArgumentParser
from pdb import set_trace as brk


def convert_field(s):
	try:
		return {"true": True, "false": False}[s]
	except KeyError:
		try:
			return int(s)
		except ValueError:
			return s


def parse_results(fname):
	with open(fname) as fp:
		line = next(iter(fp))
		line = line.strip()
		result_type = namedtuple("Result", line)

		for line in fp:
			line = line.strip()
			fields = [convert_field(f) for f in line.split()]
			yield result_type(*fields)


def parse_all_results(fname):
	ret = defaultdict(list)
	for res in parse_results(fname):
		ret[res.name].append(res)
	return ret

def score(results, field):
	x, y = [], []
	for result in results:
		x.append(result.acceptable)
		y.append(getattr(result, field))

	pb = pointbiserialr(x, y)
	return pb.correlation, pb.pvalue


def silent_counts(results):
	for k, v in results.items():
		for field in ("muts_in_sites", "total_sites", "total_singles"):
			print(k, field, score(v, field))


def make_graph_files(results):
	for k, v in results.items():
		with open("{}_true.dat".format(k), "wt") as true_fp:
			with open("{}_false.dat".format(k), "wt") as false_fp:
				for result in v:
					fp = true_fp if result.acceptable else false_fp
					print(result.count, result.max_length, file=fp)

def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	ap.add_argument("-s", "--silent-counts", action="store_true")
	ap.add_argument("-g", "--graph", action="store_true")

	args = ap.parse_args()
	results = parse_all_results(args.fname[0])

	if args.silent_counts:
		silent_counts(results)

	if args.graph:
		make_graph_files(results)

if __name__ == "__main__":
	main()

