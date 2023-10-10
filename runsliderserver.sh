#!/bin/bash
__usage="
Usage: $(basename $0) [PORT] [IPADDR]\n
\n
Options:\n
  PORT\t      Port in which to run Pluto Slider Server.\n
  IPADDR\t    IP address of server\n
"

if [[ ( $@ == "--help" ) || ( $@ == "-h" ) ]]
then
	echo -e ${__usage}
	exit 0
fi

command="
import PlutoSliderServer;

PlutoSliderServer.run_directory(\"notebooks/\", SliderServer_port=${1},SliderServer_host=\"${2}\")
"
echo ${command}
/home/gregfa/softwares/julia-1.8.5/bin/julia --optimize=0 -e "${command}"
