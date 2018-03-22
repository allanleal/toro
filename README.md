Reaktoro
========

First you need to build Reaktoro which depends on boost and cmake

module load cmake boost

In Reaktoro create a "build" directory go into it and run cmake.
If you are on Mac OS, you will need to pass the -Wno-dev flag and tell cmake where
BOOST is located.


cd build
cmake .. -Wno-dev -DBOOST_INCLUDE_DIR=/opt/moose/boost_1_63_0/include
make -j8

You should have a Reactoro library!

Toro
====

Now to build the toro MOOSE App, you need to tell MOOSE where the Reaktoro library is located.
If it's in the standard location, this line should work.

export REAKTORO_BUILD_DIR=$HOME/projects/Reaktoro/build
make -j8
