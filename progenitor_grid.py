import glob
import os
import shutil

import mesa_reader as mesa
import numpy as np
import pandas as pd


class ProgenitorGrid():
    """
    docstring
    """

    number_of_minus_models = 5
    number_of_plus_models = 3

    def __init__(self, grid_dir, output_file, output_dir):
        self.grid_dir = grid_dir
        self.output_file = output_file
        self.output_dir = output_dir
        self.log_dirs = self.find_all_log_dirs()
        self.number_of_tracks = len(self.log_dirs) * \
            (1 + self.number_of_minus_models + self.number_of_plus_models)
        self.grid = np.zeros(self.number_of_tracks, dtype=[('m_i', 'float64'),
                                                           ('rot', 'float64'),
                                                           ('z', 'float64'),
                                                           ('y', 'float64'),
                                                           ('fh', 'float64'),
                                                           ('fhe', 'float64'),
                                                           ('fsh', 'float64'),
                                                           ('mlt', 'float64'),
                                                           ('sc', 'float64'),
                                                           ('reimers', 'float64'),
                                                           ('blocker', 'float64'),
                                                           ('turbulence',
                                                            'float64'),
                                                           ('m', 'float64'),
                                                           ('model_number',
                                                            'float64'),
                                                           ('log_Teff', 'float64'),
                                                           ('log_L', 'float64'),
                                                           ('age', 'float64'),
                                                           ('m_core', 'float64'),
                                                           ('level', 'float64')])

    def find_all_log_dirs(self):
        """Returns list of LOG directories in the grid directory.

        Returns
        ----------
        list
            List of LOG directories.
        """
        return glob.glob(os.path.join(self.grid_dir, 'logs_*'))

    def evaluate_initial_grid(self):
        """Reads history files in a directory tree and creates a grid of parameters.
        Saves the grid to a file and creates a directory with .mod files.

        Parameters
        ----------

        Returns
        ----------
        """

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
            print(f"Output directory created: {self.output_dir}")
        else:
            print(f"Output directory '{self.output_dir}' already exists.")

        i = 0
        for log_dir in sorted(self.log_dirs):
            print(log_dir)
            initial_parameters = self.read_progenitor_name(
                os.path.basename(log_dir))
            data = mesa.MesaData(os.path.join(log_dir, "history.data"))
            if data.center_h1[-1] < 1e-4:
                rgb_tip = self.find_rgb_tip(data)  # model number, not index!
                default_model = rgb_tip//10*10 + 10

                for n in range(-self.number_of_minus_models, self.number_of_plus_models+1):
                    model_number = default_model + 10*n
                    self.add_one_row(i, model_number, n,
                                     initial_parameters, data)
                    shutil.copyfile(
                        os.path.join(log_dir, f"model_{model_number:05d}.mod"),
                        os.path.join(self.output_dir, self.create_mod_name(os.path.basename(log_dir), n, model_number)))
                    i += 1
            else:
                try:
                    shutil.rmtree(log_dir)
                    print(f"Deleted: {log_dir}")
                except OSError as e:
                    print(f"Error: {log_dir} : {e.strerror}")
            print('')

            self.save_grid_to_file(self.output_file)

    def save_grid_to_file(self, output):
        """Saves structured array to a text file.

        Parameters
        ----------
        output : str
            The name of output file.

        Returns
        ----------
        None
        """

        format = ''
        head = ''
        for ind, name in enumerate(self.grid.dtype.names):
            if self.grid[name].dtype == np.int64:
                format = format + '%18d '
                width = 18
            elif self.grid[name].dtype == np.float64:
                format = format + '%18.8f '
                width = 18
            elif self.grid[name].dtype == 'S80':
                format = format + '%80s '
                width = 80
            if len(name) <= width:
                if ind == 0:
                    head = head + (width - len(name)) * ' ' + name
                else:
                    head = head + (width - len(name) + 3) * ' ' + name
            else:
                print(
                    f"Column name {name} too long ({len(name)}) for a width: {width}")
                exit

        df = pd.DataFrame(self.grid)
        with open(output, 'w') as f:
            f.write(head + '\n')
        df.to_csv(output, header=False, index=False,
                  sep='\t', float_format='%18.8f', mode='a')

    def add_one_row(self, i, nr, level, initial_parameters, data):
        """Populates a single row of the grid.

        Parameters
        ----------
        i : int
            Index of the row.
        nr : int
            Model number of a progenitor.
        level : init
            Level of the model.
        initial_parameters : dict
            Initial parameters of progenitor in dict format.
        data : MesaData
            Evolutionary track (MESA history file) as MesaData object.

        Returns
        ----------
        """

        self.grid['m_i'][i] = initial_parameters['m']
        self.grid['rot'][i] = initial_parameters['rot']
        self.grid['z'][i] = initial_parameters['z']
        self.grid['y'][i] = initial_parameters['y']
        self.grid['fh'][i] = initial_parameters['fh']
        self.grid['fhe'][i] = initial_parameters['fhe']
        self.grid['fsh'][i] = initial_parameters['fsh']
        self.grid['mlt'][i] = initial_parameters['mlt']
        self.grid['sc'][i] = initial_parameters['sc']
        self.grid['reimers'][i] = initial_parameters['reimers']
        self.grid['blocker'][i] = initial_parameters['blocker']
        self.grid['turbulence'][i] = initial_parameters['turbulence']
        self.grid['m'][i] = data.star_mass[nr-1]
        self.grid['model_number'][i] = data.model_number[nr-1]
        self.grid['log_Teff'][i] = data.log_Teff[nr-1]
        self.grid['log_L'][i] = data.log_L[nr-1]
        self.grid['age'][i] = data.star_age[nr-1]
        self.grid['m_core'][i] = data.he_core_mass[nr-1]
        self.grid['level'][i] = level

    @staticmethod
    def read_progenitor_name(f_name):
        """Recovers initial values from the name of history file that
        follows specific format of progenitors' grid.

        Parameters
        ----------
        f_name : str
            Name of history file.

        Returns
        ----------
        values : dict
            Initial parameters in dict format.
        """

        s = f_name.split('_')
        values = {}
        values['m'] = np.float(s[1][1:])
        values['rot'] = np.float(s[2][3:])
        values['z'] = np.float(s[3][1:])
        values['y'] = np.float(s[4][1:])
        values['fh'] = np.float(s[5][2:])
        values['fhe'] = np.float(s[6][3:])
        values['fsh'] = np.float(s[7][3:])
        values['mlt'] = np.float(s[8][3:])
        values['sc'] = np.float(s[9][2:])
        values['reimers'] = np.float(s[10][7:])
        values['blocker'] = np.float(s[11][7:])
        values['turbulence'] = np.float(s[12][10:])
        return values

    @staticmethod
    def find_rgb_tip(data, log_Teff_max=3.8):
        """Returns model number of RGB tip.

        Parameters
        ----------
        data : MesaData
            Evolutionary track (MESA history file) as MesaData object.
        log_Teff_max : float, optional
            Maximum value of log_Teff, for which the function looks for RGB tip (max(log_L)).
            Default is 3.8.

        Returns
        ----------
        int
            Model number of RGB tip.
        """

        condition = np.logical_and(
            data.center_h1 < 1e-4, data.log_Teff < log_Teff_max)
        return data.model_number[condition][np.argmax(data.log_L[condition])]

    @staticmethod
    def create_mod_name(log_dir, level, model_number):
        """Returns a name of .mod file from directory name.

        Parameters
        ----------
        log_dir : str
            Basename of the directory name.
        level : int
            Level of the model.
        model_number : int
            Number of the model.

        Returns
        ----------
        str
            Name of progenitor's .mod file.
        """

        return f"rgb_{log_dir[5:]}_lvl{level}_{model_number}.mod"
