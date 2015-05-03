#include <stdio.h>
#include <stdlib.h>

int help_required(char *string)
{
	if (strings_equal(string, "help")) return(1);
	if (strings_equal(string, "HELP")) return(1);
	if (strings_equal(string, "-help")) return(1);
	if (strings_equal(string, "-HELP")) return (1);
	if (strings_equal(string, "--help")) return (1);
	if (strings_equal(string, "--HELP")) return (1);
	if (strings_equal(string, "-h")) return (1);
	if (strings_equal(string, "-H")) return (1);
	if (strings_equal(string, "--h")) return (1);
	if (strings_equal(string, "--H")) return (1);
	return(0);
}

void profits_to_dmdata_help()
{
	puts("------------------------------------------");
	puts("	This function to process fits file to   ");
	puts(" .dat file. Use dm and dm step to process ");
	puts("  the data and output the data as .dat file");
	puts("  -h help    To print this ");
	puts("  -o output  to name the output file if not use default");
	puts("  -n nthread default nthread = 4");
	puts("  -d dedispersion use dm to process the data");
	puts("  -z zerodm  use zerodm algorithm to process data before using dm ");
	puts("  -i dmmin   the minimum of the dm value");
	puts("  -a dmmax   the maximum of the dm value");
	puts("  -s dmstep  the dmstep of the dm value");
	puts("  -e zerodm  apply the zerodm");
	puts("------------------------------------------");
}

void get_fits_data_help()
{
	puts("-------------------------------------------");
	puts(" This program just to get raw data from the");
	puts(" psrfits files and to put all the rows to  ");
	puts(" one row, and one channel one file         ");
	puts(" the output file is ***.ch file			 ");
}

void pro_dmdata_help()
{
	puts("-------------------------------------------");
	puts(" This program to get the .dm file and then ");
	puts(" use baseline the data and output the file as .smd");
	puts(" -h help print this");
	puts("-------------------------------------------");
}

void draw_pic_help()
{
	puts("--------------------------------------------");
	puts(" This program to draw picture for float_smooth data");
	puts(" or gauss smooth data or draw peak flux on different dm in one panel");
	puts(" or draw histogram picture and so on");
	puts(" -a after the raw data minus the average of the data ");
	puts(" -r draw the raw data picture");
	puts(" -h print this ");
	puts(" -f draw float smooth data ");
	puts(" -g draw gauss smooth data ");
	puts(" -y draw gray scale map");
	puts(" -n set how many point to plot on the pgplot");
	puts("--------------------------------------------");
}

void find_peak_help()
{
	puts("-------------------------------------");
	puts(" This program only to find the pulse and");
	puts(" and draw the gray map of the pulse");
	puts(" -r read the file header");
	puts(" -d input the realdm data and find the peak sample");
	puts(" -s input the downsample parameter ");
	puts(" -h set the threshold parameter");
}
