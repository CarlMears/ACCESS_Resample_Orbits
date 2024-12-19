from interfaces import AMSR2_Lib

def main():
    for iorbit in range(680,1001):
        channel_list = [i for i in range(5,13)]
        write_swath_tbs = False
        write_swath_geoloc = False
        overwrite_existing_files = False
        #pathname_tb_output = "A/path/set/by/python/A/path/set/by/python/A/path/set/by/python/A/path/set/by/python"
        pathname_tb_output = 'L:/access/amsr2_tb_orbits/'
        AMSR2_lib_interface = AMSR2_Lib()
        AMSR2_lib_interface.process_and_resample_amsr2_l2b_orbit(
            iorbit, channel_list, write_swath_tbs, write_swath_geoloc, overwrite_existing_files, pathname_tb_output
            )
  
if __name__ == "__main__":
    main()
