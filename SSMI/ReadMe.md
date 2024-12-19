This document briefly describes the processing of SSMI V8 orbits for the ACCESS project

SSMI processing differs from the AMSR2 processing because we are starting with SSMI L1C files,
where the geolocation and conversion to brightness temperature has already been done.


So the steps are:
	Initialize the resampling (read in resampling weights, etc)
	
	Loop over orbits:
		read in L1C orbit.
		find locations of extra locations
		resample to circular footprints
		interpolate resampled Tbs to desired locations
		write out resampled orbit
		
		
For SSMI, we add 3 extra FOVs between each native low-res observations
		  we add 3 extra scans between each native low res scan