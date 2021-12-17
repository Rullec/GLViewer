#include "sim_kinect/SimKinect.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(sim_kinect, m)
{
    py::class_<cSimKinect>(m, "sim_kinect")
        .def(py::init<>())
        .def("Init", &cSimKinect::Init)
        .def("ApplyKinectHoleQuantizationNoise", &cSimKinect::ApplyKinectHoleQuantizationNoise)
        .def("SetFocalLength", &cSimKinect::SetFocalLength)
        .def("SetBaseline", &cSimKinect::SetBaseline)
        .def("GetFocalLength", &cSimKinect::GetFocalLength)
        .def("GetBaseline", &cSimKinect::GetBaseline);
}