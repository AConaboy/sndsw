# Example scripts
A set of example script to calibrate the FE and run data taking. They will not run out of the box, you need at least to adapt the configuration.
All these scripts (except `plot_dcr.py`) accept at least these two arguments:
* `--conf` or `-c` to specify a different configuration folder from `conf_example`
* `--board-ids` or `-b` to specify a subset of the DAQ boards to be used

These scripts are "tuned" to SciFi modules and should be modified to be used with different detectors.

The available scripts are:
* `disable_leds.py` enables or disables the LEDs on the board. Simple script, useful for quick checks.
* `get_status.py` another simple script, just prints the status of all the boards.
* `check_boards.py` checks that all the boards in the configuration respond. Prints the ID or names (`--names` flag) of those who don't. If all is good, a newline is printed.
* `calibration.py` used to calibrate before data taking.
* `plot_dcr.py` shows the threshold scan plots for a DAQ board.
* `set_thresholds.py` used to select the thresholds for all channels.
* `take_data.py` run the DAQ for a fixed amount of time.

These scripts are more useful to be called from a separate DAQ system:
* `start_daq.py` sets up and starts DAQ. Will fail if the DAQ is running leaving the system in a good state (i.e. taking data).
* `init_daq_extreset.py` initializes DAQ when an external reset is used. After calling this script, the reset must be sent trhough the TTCvi and then `start_daq_extreset.py` called. Will fail if the DAQ is running leaving the system in a good state (i.e. taking data).
* `start_daq_extreset.py` starts DAQ when an external reset is used. Will fail if the DAQ is running leaving the system in a good state (i.e. taking data).
* `send_heartbeat.py` must be called once per second so the event builder can build events.
* `status_daq.py` prints the status of the DAQ server and boards. To be used for monitoring. Feel free to change the printout.
* `stop_daq.py` Stops the DAQ. Will fail if the DAQ is stopped, leaving the system in a good state (i.e. not taking data).

Some scripts will disable the LEDs before running, you can disable this behavious by commenting the relevant line in each script.

## How to use these scripts
First of all, change the configuration to represent your setup.
You will definitely need to change the addresses of the boards and DAQ server (the ones set to `localhost` can stay if you run everything on the same PC).
You need to add or remove boards in the `configuration.toml` file and the corresponding folders for each board.

Once the configuration is ok, try that everything works with the `disable_leds.py` or the `get_status.py` scripts. They should give no errors. The second should print a dictionary with the status of each board.

Before taking data you need to run a full calibration, using `calibration.py`.
At the beginning you will be asked to set the bias voltage below breakdown (around 20V below is ok), do it and press enter.
Then you wll be asked to set it to the operation value to take the threshold scan.

Once this is done, use `plot_dcr.py` to check the threhold scan plots.
Verify that they make sense and use them to decide how to set the threholds.
Thresholds are then set with `set_thresholds.py t1 t2` (by default it sets T1 to a fixed value, T2 to a fixed rate and E to 62. You can modify this in the script).

At this point you are ready to take data.
Use the `take_data.py <run> [-t <duration>]` script. `<run>` is the desired run number (set it to -1 to have it automatically increased by one from the last run) and `<duration>` is how many seconds the DAQ will last.
The script can be modified to have different type of requirements, like a fixed number of events per run and so on.