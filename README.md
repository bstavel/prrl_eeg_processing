# PRRL EEG Anlaysis

These scripts are used to do eeg preprocessing for the PRRL study. The preprocessing is optimized for
ERP analysis and are based off scripts from Anne Collins, Sarah Master, and Amy Zou.

## Acquistion Assumptions

After informed consent procedures and demographic data collection, each participant should be situated in our BioSemi EEG equipment. To begin achieve this, subjects will be seated in a comfortable chair in the preparation room. Single-use alcohol preparation pads will be used to remove sweat and dried skin from seven locations in which non-capelectrodes are to be placed.  These locations will include: the outer canthi of the rightand left eyes, above and below the center of the right eye (to track electro-oculographicactivity), on the left and right mastoids, and on the bridge of the nose toward the nosetip. The participant’s head circumference and inion–nasion distance will be then measuredand, based on the circumference, an appropriately-sized elastic cap will be selected. Afterplacing any non-cap electrodes expected to lie underneath the cap (e.g., on the mastoidsand above the right eye) using small adhesive disks, the cap will be placed according tothe 10–20 electrode location system. Next, conductive electrode gel will be placed in the 64 electrode sockets in the cap. Then the 64 scalp electrodes will be connected to the capand, finally, any remaining non-cap electrodes will be put in place

_Documentation of these steps was provided by Beth Baribault_

## Basic Preprocessing Steps

0. Load the data into the EEGLAB software
1. Full walk through of the EEG timeseries, removing epochs manually via the EEGLAB GUI.
2. Option to interpolate bad channels.
3. ICA and removal of ICA components that contain aritifacts such as eye blinks
4. Filtering between 0.5-20.
5. Visualization of epochs by channel, resulting in 64 figures (one for each electrode), and option to more closely inspect the epochs in a given channel.
6. Displays the most divergent epochs for channels selected for closer inspection
7. Option to visually inspect divergent epochs and mark them for removal.
8. Option to mark additional bad channels.
9. Any additional bad epochs or channels are fed back into steps 1 and 2.
10. Steps 3-8 are then repreated.

The current script is built to accommodate no more than 3 passes through the data without manual manipulations of the script.

## PRRL Specific Information

EEG Flag descriptions can be found here: https://docs.google.com/spreadsheets/d/19s4jlUIBumjqCz5WKiB08tHt3FIf2Kv3-DGEBlDefS4/edit?usp=sharing


