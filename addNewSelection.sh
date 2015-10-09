#!/bin/bash

if [ ! -e cutflow${1}.cpp ]; then 
  echo "#include \"LLGAnalysis.h\"" >> cutflow${1}.cpp
  echo "void LLGAnalysis::Setup${1}() {" >> cutflow${1}.cpp
  echo  >> cutflow${1}.cpp
  echo "    // setup the cutflow">> cutflow${1}.cpp
  echo  >> cutflow${1}.cpp
  echo "    // and the histograms" >> cutflow${1}.cpp
  echo  >> cutflow${1}.cpp
  echo "    return;">> cutflow${1}.cpp
  echo "}">> cutflow${1}.cpp
  echo  >> cutflow${1}.cpp
  echo "void LLGAnalysis::${1}Selection() {">> cutflow${1}.cpp
  echo  >> cutflow${1}.cpp
  echo "    return;">> cutflow${1}.cpp
  echo "}">> cutflow${1}.cpp
else
  echo "SOURCE FILE ALREADY EXISTS. NOT CREATING A NEW SOURCE FILE FOR Region ${1}"
fi




if ! grep -q "Setup${1}()" LLGAnalysis.h; then
  sed -ie '/INSERT YOUR SELECTION HERE/ i\        void Setup'${1}'();' LLGAnalysis.h
  sed -ie '/INSERT YOUR SELECTION HERE/ i\        void '${1}'Selection();' LLGAnalysis.h
else
  echo "Setup${1} ALREADY KNOW TO LLGAnalysis.h. LEAVING FILE UNCHANGED"
fi

if ! grep -q "Setup${1}()" LLGAnalysis.cpp; then
  sed -ie '/SETUP YOUR SELECTION HERE/ i\    else if( SELECTION == "'${1}'" ) Setup'${1}'();' LLGAnalysis.cpp
  sed -ie '/CALL YOUR SELECTION HERE/ i\        else if( SELECTION == "'${1}'" ) '${1}'Selection();' LLGAnalysis.cpp
else
  echo "Setup${1} ALREADY KNOW TO LLGAnalysis.cpp. LEAVING FILE UNCHANGED"
fi


if ! grep -q "cutflow${1}.o" Makefile; then 
  sed -ie 's/cutflowSignalRegion.o/cutflowSignalRegion.o cutflow'${1}'.o/' Makefile
else
  echo "cutflow${1}.o ALREADY KNOWN TO Makefile. LEAVING FILE UNCHANGED"
fi
