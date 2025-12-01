import subprocess
import csv
import io
import matplotlib.pyplot as plt

EXE_PATH = r'out\build\mainpreset\Debug\SimplePKPDProject.exe'
CSV_PATH = r"pk_output.csv"
ARGS = [
    EXE_PATH,
    "--dose", "100",
    "--t12", "13.86",
    "--f", "1",
    "--tau", "24",
    "--ndoses", "6",
    "--ka", "0.10",
]

def run_simulation():
    """
    """
    result = subprocess.run(
        ARGS,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    print(result.stdout)
    if result.stderr:
        print(result.stderr)

def parse_csv(csv_path):
    """
    Reads the CSV file produced by the exe.
    Adjust column names as needed.
    """
    times = []
    concs = []

    with open(csv_path, "r", encoding="utf-8") as f:
        first = f.readline().strip()

        if first.lower().startswith("sep="):
            delimiter = first.split("=", 1)[1]  # usually ";"
            # Next line will be the header
        else:
            delimiter = ";"
            f.seek(0)  # rewind if no sep= line

        # reader = csv.DictReader(f)
        reader = csv.DictReader(f, delimiter=delimiter)

        for row in reader:
            # Convert decimal commas → dots
            t = float(row["time"].replace(",", "."))
            c = float(row["Ac"].replace(",", "."))   # adjust column if needed

            times.append(t)
            concs.append(c)

    return times, concs

# def parse_csv(csv_text):
#     """
#     Parses CSV text into lists of time and concentration.
#     Assumes columns are named 'time' and 'conc' in the header.
#     """
#     times = []
#     concs = []

#     f = io.StringIO(csv_text)
#     print(f"f={f}")
#     reader = csv.DictReader(f)
#     for row in reader:
#         print(f"row: {row}")
#         # adjust names to match your program's header
#         t = float(row["time"])
#         c = float(row["Ac"])
#         times.append(t)
#         concs.append(c)

#     return times, concs

def plot_concentration_time(times, concs):
    plt.figure()
    plt.plot(times, concs) # marker="o"
    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.title("Concentration vs Time")
    plt.grid(True)
    plt.show()   # non-blocking
    # plt.pause(0.001)        


def main():
    run_simulation()
    times, concs = parse_csv(CSV_PATH)
    
    for c in concs:
        print(c)
    
    plot_concentration_time(times, concs)

    # input("Press Enter to exit...")

if __name__ == "__main__":
    main()
