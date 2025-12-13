import subprocess
import csv
import io
import matplotlib.pyplot as plt

EXE_PATH = r'out\build\mainpreset\Debug\SimplePKPDProject.exe'
CSV_PATH = r"pk_output.csv"
ARGS = [
    EXE_PATH, 
    # "--out", # CSV_PATH,
    "--dose", "100",
    "--t12", "12",
    # "--f", "1",
    # "--tau", "24",
    # "--ndoses", "2",
    # "--ka", "0.15",
    # "--ev"
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
        print("Debug msg:")
        print(result.stderr)
        print("-----")

    return result

def parse_csv(csv_path):
    """
    Reads the CSV file produced by the exe.
    Adjust column names as needed.
    """
    times   = []
    concs   = []
    depot_c = []

    with open(csv_path, "r", encoding="utf-8") as f:
        first = f.readline().strip()

        if first.lower().startswith("sep="):
            delimiter = first.split("=", 1)[1]  # usually ";"
            # Next line will be the header
        else:
            delimiter = ";"
            f.seek(0)  # rewind if no sep= line

        reader = csv.DictReader(f, delimiter=delimiter)

        for row in reader:
            # Convert decimal commas → dots
            t = float(row["time"].replace(",", "."))
            d = float(row["Ag"].replace(",", "."))
            c = float(row["Ac"].replace(",", "."))

            times.append(t)
            depot_c.append(d)
            concs.append(c)

    return times, concs, depot_c

def read_csv_cout(csv_text):
    f = io.StringIO(csv_text)
    first = f.readline().strip()
    
    if first.lower().startswith("sep="):
        delimiter = first.split("=", 1)[1]
    else:
        delimiter = ";" 
        f.seek(0)

    reader = csv.DictReader(f, delimiter=delimiter, skipinitialspace=True)
    
    if reader.fieldnames:
        reader.fieldnames = [h.strip().lstrip("\ufeff") for h in reader.fieldnames]
    has_ag = reader.fieldnames and ("Ag" in reader.fieldnames)

    times, concs= [], []
    depot_c = [] if has_ag else None

    for row in reader:
        t = float(row["time"].strip().replace(",", "."))
        c = float(row["Ac"].strip().replace(",", "."))
        if has_ag:
            depot_c.append(float(row["Ag"].strip().replace(",", ".")))

        times.append(t); concs.append(c)

    return times, concs, depot_c

def plot_concentration_time(times, concs, depot_c=None):
    plt.figure()
    plt.plot(times, concs, label="Central")
    if (depot_c):
        plt.plot(times, depot_c, label="Depot")
    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.title("Concentration vs Time")
    plt.grid(True)
    plt.legend()
    plt.show()

def main():
    result = run_simulation()

    if ("-o" in ARGS) or ("--out" in ARGS): 
        print("Reading .csv file output")
        times, concs, depot_c = parse_csv(CSV_PATH)
    else: 
        print("Reading stdcout")
        times, concs, depot_c = read_csv_cout(result.stdout)

    plot_concentration_time(times, concs, depot_c)
    

    # input("Press Enter to exit...")

if __name__ == "__main__":
    main()
