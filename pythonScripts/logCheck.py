from pathlib import Path

directory = "/Users/mcgurn/Downloads/tmp_13698594"

files = Path(directory).glob('*')
for file in files:
    with open(file) as openFile:

        line = openFile.readlines()[-18]
        if "Event begin: Radiation::RadReturn" not in line:
            print(file)
            print(line)
