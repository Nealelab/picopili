#! /bin/sh

###########
#
# get_refs.sh
# Retrives reference files
#
# - ricopili platform guessing (plague) files
# - ricopili build guessing (buigue) files
#
###########

echo " "
echo "### External reference file downloader for picopili ###"
echo " "
echo "Picopili depends on a few curated reference files"
echo "from ricopili. If ricopili is installed on your"
echo "platform, will set up symbolic links to the required"
echo "files. Otherwise, will download the files."
echo " "

echo "### BEGIN ###"
echo " "

# setup
rp_conf="$HOME/ricopili.conf"
SERVER="https://personal.broadinstitute.org/rwalters/picopili_files"
SCRIPT=$(readlink -f "$0")
BINLOC=$(dirname "$SCRIPT")
LIBLOC=`echo $(dirname "$BINLOC")"/lib"`
rp=false
li_done=false
hm_done=false

if [ -d "$LIBLOC/buigue" ]; then
	echo "WARNING: Found existing folder $LIBLOC/buigue. Contents may be overwritten."
	echo "(pausing to allow cancel...)"
	sleep 3
	echo "(continuing)"
else
	mkdir "$LIBLOC/buigue"
fi

if [ -d "$LIBLOC/plague" ]; then
	echo "WARNING: Found existing folder $LIBLOC/plague. Contents may be overwritten."
	echo "(pausing to allow cancel...)"
	sleep 3
	echo "(continuing)"
else
	mkdir "$LIBLOC/plague"
fi


# check/read config file
if [ -e "$rp_conf" ]; then
	rp=true
	echo "Found existing ricopili configuration. Reading..."
	liloc=`awk '$1=="liloc"{print $2}' $rp_conf`
	hmloc=`awk '$1=="hmloc"{print $2}' $rp_conf`
else
	echo "No ricopili configuration found."
fi

lifiles=("snp.txt.pos.scz49.gz" "snp125.txt.pos.scz49.gz" "snp130.txt.pos.scz49.gz" "snp138.txt.pos.scz49.gz" "last")
hmfiles=("snp_platform_collection.txt.new.0815.gz" "snp_platform_collection.txt.new.0416a.gz" "snp_platform_collection.txt.new.0114.gz" "last")

# link creation from ricopili references
if [ "$rp" = 'true' ]; then 

	if [ -d "$liloc" ]; then
		
		for finame in ${lifiles[@]}; do
			
			if [ "$finame" = "last" ]; then
				li_done=true
			else
				echo "$liloc/$finame"
				ln -sfn "$liloc/$finame" "$LIBLOC/buigue" || break
			fi
		done
	fi
	
	if [ "$li_done" = 'false' ]; then
		echo "Failed to link all files from liftover directory $liloc"
	fi


	if [ -d "$hmloc" ]; then
		for finame in ${hmfiles[@]}; do
			
			if [ $finame == "last" ]; then
				hm_done=true
			else
				echo "$hmloc/$finame"
				ln -sfn "$hmloc/$finame" "$LIBLOC/plague" || break
			fi
		done
	fi

	if [ "$hm_done" = 'false' ]; then
		echo "Failed to link all platform references from directory $hmloc"
	fi
fi

if [ "$li_done" = 'false' ]; then 
	to_dl=true
elif [ "$hm_done" = 'false' ]; then
	to_dl=true
else
	to_dl=false
fi

# wget external
if [ "$to_dl" = 'true' ]; then

	# warn of internet access
	echo " "
	echo "WARNING: Preparing to download reference files from:"
	echo "$SERVER"
	echo " "
	echo "Expected total file size is ~300 MB, minus existing"
	echo "files already linked/downloaded."
	echo " "
	echo "If you do not have web access, or if you do not want"
	echo "to download these files now, please cancel now."
	echo " "
	echo "Will begin in 10 sec..."
	echo " "
	sleep 10

	for finame in ${lifiles[@]}; do
		
		if [ "$finame" = "last" ]; then
			continue
		else
			echo " "
			echo "Next file: $SERVER/$finame"
			# wget --no-check-certificate "$SERVER/$finame" "$LIBLOC/buigue/$finame"
			curl -o "$LIBLOC/buigue/$finame" "$SERVER/$finame"
		fi
	done
	for finame in ${hmfiles[@]}; do

		if [ "$finame" = "last" ]; then
			continue
		else
			echo " "
			echo "Next file: $SERVER/$finame"
			# wget --no-check-certificate "$SERVER/$finame" "$LIBLOC/plague/$finame"
			curl -o "$LIBLOC/plague/$finame" "$SERVER/$finame"
		fi
	done
fi

echo " "
echo "### Finished ###"
echo " "

# eof
