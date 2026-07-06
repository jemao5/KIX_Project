import glob
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_glob", help="glob to find pckl files to be merged")
    parser.add_argument("output_path", help="path for merged pckl file")

    return parser.parse_args()

def main():
    args = parse_args()

    # find all per-chunk pickles
    chunk_files = sorted(glob.glob(args.input_glob))
    print(f"Found {len(chunk_files)} chunk pickles")

    all_data = []
    for cf in chunk_files:
        with open(cf, "rb") as f:
            chunk_data = pickle.load(f)   # a list of interaction dicts
        all_data.extend(chunk_data)        # append this chunk's entries to the master list
        print(f"  {cf}: {len(chunk_data)} entries")

    print(f"Total merged entries: {len(all_data)}")

    # write the merged pickle
    with open(args.output_path, "wb") as f:
        pickle.dump(all_data, f)
    print("Wrote merged pkl")


if __name__ == "__main__":
    main()
