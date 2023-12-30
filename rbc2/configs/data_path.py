from pathlib import Path

path_to_data_folder = str(Path(__file__).parents[1]) + '/data'

if __name__ == '__main__':
    print(path_to_data_folder)

