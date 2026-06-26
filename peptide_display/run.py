from load_raw_data import load_raw_data
from generate_explore_indexes import generate_explore_indexes
from generate_website_indexes import generate_website_indexes


def main():
    print("Hello from peptide-display part of the pipeline!\n")
    load_raw_data()
    generate_explore_indexes()
    generate_website_indexes()


if __name__ == "__main__":
    main()
