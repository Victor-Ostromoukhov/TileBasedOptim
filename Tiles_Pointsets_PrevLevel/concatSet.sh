#!/usr/bin/env zsh

#for n in  03 09 27 81 243 729 2187 6561 19683 59049

foreach i  (*/*03.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_3.dat
  echo "#" >> pts_3.dat
end

foreach i  (*/*09.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_9.dat
  echo "#" >> pts_9.dat
end

foreach i  (*/*27.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_27.dat
  echo "#" >> pts_27.dat
end

foreach i  (*/*81.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_81.dat
  echo "#" >> pts_81.dat
end

foreach i  (*/*243.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_243.dat
  echo "#" >> pts_243.dat
end

foreach i  (*/*729.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_729.dat
  echo "#" >> pts_729.dat
end

foreach i  (*/*2187.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_2187.dat
  echo "#" >> pts_2187.dat
end


foreach i  (*/*6561.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_6561.dat
  echo "#" >> pts_6561.dat
end


foreach i  (*/*19683.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_19683.dat
  echo "#" >> pts_19683.dat
end


foreach i  (*/*59049.dat)
  echo $i
  cat $i | awk '{print $2" "$3}' >> pts_59049.dat
  echo "#" >> pts_59049.dat
end
