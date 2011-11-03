%module lsd
%module mymodule 
%{
#include "lsd.h"
%}

%include "lsd.h"

%include "carrays.i"
%array_class(double, doubleArray);

