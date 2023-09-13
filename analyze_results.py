from scipy.stats import pointbiserialr
from collections import namedtuple, defaultdict
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
	with open("test.txt") as fp:
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


def main():
	results = parse_all_results("test.txt")

	for k, v in results.items():
		for field in ("muts_in_sites", "total_sites", "total_singles"):
			print(k, field, score(v, field))


if __name__ == "__main__":
	main()

