#
# after generating the verification pages ("reports") symlink them into a
# group-specific subdirectory.  errors in the symlinking (e.g., if the symlink
# already exists) are sent to /dev/null.
#
targetdirs=`ls -d ????????????????*/`

for g in $targetdirs; do

  echo 'symlinking '$g' ...'

  outdir=/Users/luke/Dropbox/proj/cdips_followup/results/followup_planning/2020A_youngstar

  tolink=`ls $g*png`

  # trim off the last trailing "/"
  prestr=`echo $g | sed 's/.$//'`

  for t in $tolink; do
    echo $t
    
    str0=`dirname $t`
    str1=`basename $t`

    ln -s `realpath $t` "$outdir"/"$str0"_"$str1"
  done

done
