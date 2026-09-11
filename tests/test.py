import sys
from effspm import HTMiner as mine

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <data_file> <minsup>")
        sys.exit(1)

    data_file = sys.argv[1]
    try:
        minsup = float(sys.argv[2])
    except ValueError:
        print("minsup must be a number (e.g. 0.01)")
        sys.exit(1)

    # Call the C++ binding wrapper
    result = mine(data_file, minsup)
    
    # Extract patterns list and execution time from the returned dictionary
    pattern_list = result.get('patterns', [])
    exec_time = result.get('time', 0.0)

    print(f"Found {len(pattern_list)} patterns in {exec_time:.4f} seconds\n")
    
    # Print the first 10 mined patterns
    for i, pat in enumerate(pattern_list[:10]):
        print(f"Pattern {i+1}: {pat}")

if __name__ == "__main__":
    main()