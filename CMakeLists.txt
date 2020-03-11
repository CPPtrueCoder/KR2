cmake_minimum_required(VERSION 3.15)
project(KR)

set(CMAKE_CXX_STANDARD 14)

add_executable(KR main.cpp Object.cc Object.h DormandPrinceIntegrator.cc
        DormandPrinceIntegrator.h arenstorf.h arenstorf.cc)