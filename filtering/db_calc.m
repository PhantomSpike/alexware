function db_level = db_calc(data)
xrms = rms(data); %Calculate rms of the signal
Q = 2*(10^-5); %The constant to which all sound levels are compared ~ 20uPa to get the result in sound pressure level (SPL)
db_level = mag2db(xrms/Q); %Calculate the sound level in dB sound pressure level 