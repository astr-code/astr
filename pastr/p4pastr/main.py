import sys
import os

def monitor_process(filename):
    from . import _pastr
    import matplotlib.pyplot as plt

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    print("Monitoring file:", filename)

    if len(sys.argv) < 4:
        print("Usage: python -m p4pastr monitor <filename> <number>")
        sys.exit(1)

    filename = sys.argv[2]

    try:
        number = int(sys.argv[3])   # or float(sys.argv[3])
    except ValueError:
        raise ValueError("The third argument must be a number")

    ncol, nrecl = _pastr.kernels.monitor_data_count(filename)
    time,var=_pastr.kernels.read_monitor_data(filename,nrecl,ncol,number)

    plt.plot(time, var)
    plt.xlabel("Time")
    plt.ylabel("Physical variable")
    plt.grid(True)
    plt.show()

    print("Job run_pastr finished")

def run_pastr():
    print("  ** currently nothing, try:")
    print("  python -m p4pastr.main monitor <filename> <number>")

def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python -m p4pastr.main monitor <filename> <number>")
        print("  python -m p4pastr.main <other-command>")
        sys.exit(1)

    command = sys.argv[1]

    if command == "monitor":
        if len(sys.argv) < 3:
            print("Error: monitor requires a filename")
            sys.exit(1)

        filename = sys.argv[2]
        monitor_process(filename)
    else:
        run_pastr(command)

if __name__ == "__main__":
    main()
