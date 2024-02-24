from scipy.stats import pointbiserialr, kstest
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
		while True:
			line = next(iter(fp))
			line = line.strip()
			if line.startswith('#'):
				continue

			result_type = namedtuple("Result", line)
			break

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


def rates(results, max_count=None,
		  exact_count=None, require_not_interleaved=False):
	for k, v in results.items():
		good, total = 0, 0
		for result in v:
			acceptable = result.acceptable

			# The "acceptable" field in the results is redundant but it's still
			# worth working it out in the Go program just so we can print it
			# out as we go along. So we might as well check it here.
			assert acceptable == (result.unique and result.max_length < 8000)

			if max_count is not None:
				acceptable = acceptable and result.count <= max_count

			if exact_count is not None:
				acceptable = acceptable and result.count == exact_count

			if require_not_interleaved:
				acceptable = acceptable and not result.interleaved

			if acceptable:
				good += 1
			total += 1

		print("{}: {}/{} {:.4g}%".format(k, good,
			total, float(good * 100) / total))


def normalized_positions(result):
	positions = [float(x) for x in result.positions.strip('[]').split(',')]
	return [x / float(result.genome_len) for x in positions]


def ks(positions):
	return kstest(positions, "uniform").pvalue


def mean_abs_diff(positions):
	n = len(positions) + 1
	d = 0
	for i, pos in enumerate(positions):
		d += abs(pos - (i+1) / n)
	return d * 1/(n-1)


def check_results(results, test, ref_positions):
	ref = test(ref_positions)
	greater = 0
	total = 0.0

	print("Reference value: {:.3f}".format(ref))
	print("Starting point, mean value, % greater than reference")

	for k, v in results.items():
		greater = 0
		total = 0.0
		for result in v:
			positions = normalized_positions(result)
			val = test(positions)
			total += val
			if val > ref:
				greater += 1
		average = total / len(v)
		print("{} {:.3f} {:.3f}%".format(k, average, (greater * 100) / len(v)))
	print()


def added_removed(results):
	print("Average numbers of sites added and removed")
	for k, v in results.items():
		total_added, total_removed = 0, 0
		for result in v:
			total_added += result.added
			total_removed += result.removed

		n = len(v)
		added = total_added / n
		removed = total_removed / n
		print("{}: {} muts {:.2f} added {:.2f} removed".format(k,
		   v[0].num_muts, added, removed))
	print()


def main():
	ap = ArgumentParser()
	ap.add_argument("fname", nargs=1)
	ap.add_argument("-s", "--silent-counts", action="store_true")
	ap.add_argument("-g", "--graph", action="store_true")
	ap.add_argument("-r", "--rates", action="store_true", default=True)
	ap.add_argument("-m", "--max-count", type=int)
	ap.add_argument("-e", "--exact-count", type=int)
	ap.add_argument("-i", "--require-not-interleaved", action="store_true")

	args = ap.parse_args()
	results = parse_all_results(args.fname[0])

	if args.silent_counts:
		silent_counts(results)

	if args.graph:
		make_graph_files(results)

	if args.max_count or args.exact_count or args.require_not_interleaved:
		rates(results, args.max_count,
			args.exact_count, args.require_not_interleaved)
	elif args.rates:
		rates(results)

	print()
	added_removed(results)

	WH1 = [0.0733036819048256, 0.32605424204929273,
		0.57947363140822, 0.6009764906531118, 0.8059726448851285]

	print("Kolomogorov Smirnov")
	check_results(results, ks, WH1)

	print("Mean Absolute Difference from perfectly uniform")
	check_results(results, mean_abs_diff, WH1)


if __name__ == "__main__":
	main()
