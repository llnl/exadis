/*---------------------------------------------------------------------------
 *
 *	ExaDiS python binding module
 *
 *	Nicolas Bertin
 *	bertin1@llnl.gov
 *
 *-------------------------------------------------------------------------*/

#pragma once
#ifndef EXADIS_PYBIND_H
#define EXADIS_PYBIND_H

#include <iostream>
#include <exadis.h>
#include <driver.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/*---------------------------------------------------------------------------
 *
 *    Struct:        type_caster
 *
 *-------------------------------------------------------------------------*/
namespace pybind11 { namespace detail {
    
    template <> struct type_caster<ExaDiS::Vec3> {
    public:
        /*
         * This macro establishes the name 'Vec3' in
         * function signatures and declares a local variable
         * 'value' of type Vec3
         */
        PYBIND11_TYPE_CASTER(ExaDiS::Vec3, _("Vec3"));

        /*
         * Conversion part 1 (Python->C++): convert a PyObject into a Vec3
         * instance or return false upon failure. The second argument
         * indicates whether implicit conversions should be applied.
         */
        bool load(handle src, bool) {
            if (!isinstance<sequence>(src)) return false;
            sequence seq = reinterpret_borrow<sequence>(src);
            if (seq.size() != 3)
                throw value_error("Expected sequence of length 3 for Vec3 type");
            value[0] = seq[0].cast<double>();
            value[1] = seq[1].cast<double>();
            value[2] = seq[2].cast<double>();
            return true;
        }

        /*
         * Conversion part 2 (C++ -> Python): convert a Vec3 instance into
         * a Python object. The second and third arguments are used to
         * indicate the return value policy and parent object (for
         * ``return_value_policy::reference_internal``) and are generally
         * ignored by implicit casters.
         */
        static handle cast(ExaDiS::Vec3 src, return_value_policy, handle) {
            return py::make_tuple(src.x, src.y, src.z).release();
        }
    };
    
    template <> struct type_caster<ExaDiS::Vec3i> {
    public:
        PYBIND11_TYPE_CASTER(ExaDiS::Vec3i, _("Vec3i"));
        bool load(handle src, bool) {
            if (!isinstance<sequence>(src)) return false;
            sequence seq = reinterpret_borrow<sequence>(src);
            if (seq.size() != 3)
                throw value_error("Expected sequence of length 3 for Vec3i type");
            value[0] = seq[0].cast<int>();
            value[1] = seq[1].cast<int>();
            value[2] = seq[2].cast<int>();
            return true;
        }
        static handle cast(ExaDiS::Vec3i src, return_value_policy, handle) {
            return py::make_tuple(src.x, src.y, src.z).release();
        }
    };
    
    template <> struct type_caster<ExaDiS::Mat33> {
    public:
        PYBIND11_TYPE_CASTER(ExaDiS::Mat33, _("Mat33"));
        bool load(handle src, bool) {
            if (!isinstance<sequence>(src)) return false;
            sequence seq = reinterpret_borrow<sequence>(src);
            if (seq.size() == 3) {
                value[0] = seq[0].cast<ExaDiS::Vec3>();
                value[1] = seq[1].cast<ExaDiS::Vec3>();
                value[2] = seq[2].cast<ExaDiS::Vec3>();
            } else if (seq.size() == 9) {
                value[0] = ExaDiS::Vec3(seq[0].cast<double>(), seq[1].cast<double>(), seq[2].cast<double>());
                value[1] = ExaDiS::Vec3(seq[3].cast<double>(), seq[4].cast<double>(), seq[5].cast<double>());
                value[2] = ExaDiS::Vec3(seq[6].cast<double>(), seq[7].cast<double>(), seq[8].cast<double>());
            } else {
                throw value_error("Expected sequence of length 9 or 3x3 for Mat33 type");
            }
            return true;
        }
        static handle cast(ExaDiS::Mat33 src, return_value_policy, handle) {
            return py::make_tuple(src.rowx, src.rowy, src.rowz).release();
        }
    };
    
    template <> struct type_caster<ExaDiS::NodeTag> {
    public:
        PYBIND11_TYPE_CASTER(ExaDiS::NodeTag, _("NodeTag"));
        bool load(handle src, bool) {
            if (!isinstance<sequence>(src)) return false;
            sequence seq = reinterpret_borrow<sequence>(src);
            if (seq.size() != 2)
                throw value_error("Expected sequence of length 2 for NodeTag type");
            value.domain = seq[0].cast<int>();
            value.index  = seq[1].cast<int>();
            return true;
        }
        static handle cast(ExaDiS::NodeTag src, return_value_policy, handle) {
            return py::make_tuple(src.domain, src.index).release();
        }
    };
    
}} // namespace pybind11::detail

namespace ExaDiS { namespace pybind {

/*---------------------------------------------------------------------------
 *
 *    Utility functions
 *
 *-------------------------------------------------------------------------*/
void initialize(int num_threads=-1, int device_id=0, bool verbose=true);
void finalize();
std::vector<int> map_node_tags(DeviceDisNet* net, std::vector<NodeTag>& tags);
void set_positions(System* system, std::vector<Vec3>& pos);
void set_forces(System* system, std::vector<Vec3>& forces, std::vector<NodeTag>& tags);
void set_velocities(System* system, std::vector<Vec3>& vels, std::vector<NodeTag>& tags);
std::vector<Vec3> get_positions(System* system);
std::vector<Vec3> get_forces(System* system);
std::vector<Vec3> get_velocities(System* system);


/*---------------------------------------------------------------------------
 *
 *    ExaDisNet binding
 *    Wrap the dislocation network into an ExaDiS system object.
 *    This allows to save on overhead time when driving a GPU simulation
 *    and prevents unnecessary memory copies between spaces.
 *
 *-------------------------------------------------------------------------*/
struct ExaDisNet {
    System* system = nullptr;
    
    ExaDisNet() {
        system = make_system(new SerialDisNet(), Crystal(), Params());
        system->pyexadis = true;
    }
    
    ExaDisNet(System* _system) : system(_system) {
        system->pyexadis = true;
    }
    
    ExaDisNet(Cell& cell,
              std::vector<std::vector<double> >& nodes_array, 
              std::vector<std::vector<double> >& segs_array)
    {
        SerialDisNet* net = new SerialDisNet(cell);
        net->set_nodes_array(nodes_array);
        net->set_segs_array(segs_array);
        net->sanity_check();
        system = make_system(net, Crystal(), Params());
        system->pyexadis = true;
    }
    
    void import_data(Cell& cell,
                     std::vector<std::vector<double> >& nodes_array, 
                     std::vector<std::vector<double> >& segs_array)
    {
        if (!system)
            ExaDiS_fatal("Error: cannot import data in unitialized ExaDisNet object\n");
        SerialDisNet* net = system->get_serial_network();
        net->cell = cell;
        net->set_nodes_array(nodes_array);
        net->set_segs_array(segs_array);
        net->sanity_check();
        net->update();
    }
    
    System* adjust_system(Params& params) {
        system->params = params;
        if (system->crystal != params.crystal)
            system->crystal = Crystal(params.crystal);
        return system;
    }
    
    int number_of_nodes() { return system->Nnodes_total(); }
    int number_of_segs() { return system->Nsegs_total(); }
    bool is_sane() { return system->get_serial_network()->sanity_check(); }
    
    Cell get_cell() { return system->get_serial_network()->cell; }
    std::vector<std::vector<double> > get_nodes_array() { return system->get_serial_network()->get_nodes_array(); }
    std::vector<std::vector<double> > get_segs_array() { return system->get_serial_network()->get_segs_array(); }
    std::vector<Vec3> get_forces() { return ExaDiS::pybind::get_forces(system); }
    std::vector<Vec3> get_velocities() { return ExaDiS::pybind::get_velocities(system); }
    py::tuple get_plastic_strain() { return py::make_tuple(system->dEp, system->dWp, system->density); }
    
    void set_positions(std::vector<Vec3>& pos) { ExaDiS::pybind::set_positions(system, pos); }
    void set_forces(std::vector<Vec3>& forces, std::vector<NodeTag>& tags) { ExaDiS::pybind::set_forces(system, forces, tags); }
    void set_velocities(std::vector<Vec3>& vels, std::vector<NodeTag>& tags) { ExaDiS::pybind::set_velocities(system, vels, tags); }
    
    Crystal* get_crystal() { return &system->crystal; }
    SerialDisNet* get_serial_network() { return system->get_serial_network(); }
    
    SerialDisNet::DisLinks physical_links() { return system->get_serial_network()->physical_links(); }
    
    void write_data(std::string filename) { system->get_serial_network()->write_data(filename); }
};

struct SystemBind : ExaDisNet {
    SystemBind(ExaDisNet disnet, Params params) : ExaDisNet()
    {
        SerialDisNet* net = disnet.system->get_serial_network();
        system = make_system(net, Crystal(params.crystal), params);
        system->params.check_params();
    }
    void set_neighbor_cutoff(double cutoff) {
        system->neighbor_cutoff = cutoff;
    }
    void set_applied_stress(std::vector<double> applied_stress) { 
        system->extstress = Mat33().voigt(applied_stress.data()); 
    }
    void print_timers(double timetot, bool dev) { system->print_timers(timetot, dev); }
};


/*---------------------------------------------------------------------------
 *
 *    Force binding
 *
 *-------------------------------------------------------------------------*/
struct ForceBind {
    enum ForceModel {
        LINE_TENSION_MODEL, CUTOFF_MODEL, DDD_FFT_MODEL, 
        SUBCYCLING_MODEL, GLOBAL_MODEL, PYTHON_MODEL,
        FORCE_FFT,
    };
    Force* force = nullptr;
    int model = -1;
    Params params;
    double neighbor_cutoff = 0.0;
    bool pre_computed = false;
    ForceBind(Force* _force, int _model, Params _params, double cutoff=0.0) : 
    force(_force), model(_model), params(_params), neighbor_cutoff(cutoff) {}
    
    ForceFFT* get_force_fft() {
        ForceFFT* forcefft = nullptr;
        if (model == DDD_FFT_MODEL) {
            forcefft = static_cast<ForceType::DDD_FFT_MODEL*>(force)->get_force2()->get_flong();
        } else if (model == GLOBAL_MODEL) {
            forcefft = static_cast<ForceType::GLOBAL_MODEL*>(force)->get<ForceGlobal::FORCE_FFT>();
        } else if (model == FORCE_FFT) {
            forcefft = static_cast<ForceFFT*>(force);
        } else { 
            ExaDiS_fatal("Error: get_force_fft() method requires DDD_FFT_MODEL or FORCE_FFT model\n");
        }
        return forcefft;
    }
    void pre_compute(SystemBind& sysbind) { force->pre_compute(sysbind.system); }
    void compute(SystemBind& sysbind) { force->compute(sysbind.system); }
    
    void pre_compute_force(ExaDisNet& disnet) {
        System* system = disnet.adjust_system(params);
        force->pre_compute(system);
        pre_computed = true;
    }
    std::vector<Vec3> compute_force(ExaDisNet& disnet, std::vector<double> applied_stress,
                                    bool pre_compute) {
        System* system = disnet.adjust_system(params);
        system->extstress = Mat33().voigt(applied_stress.data());
        if (pre_compute) {
            force->pre_compute(system);
            pre_computed = true;
        }
        force->compute(system);
        std::vector<Vec3> forces = get_forces(system);
        return forces;
    }
    Vec3 compute_node_force(ExaDisNet& disnet, int i, 
                            std::vector<double> applied_stress) {
        System* system = disnet.adjust_system(params);
        system->extstress = Mat33().voigt(applied_stress.data());
        // Warning: the user must ensure the pre_compute is up-to-date...
        if (!pre_computed) {
            force->pre_compute(system);
            pre_computed = true;
        }
        Vec3 f = force->node_force(system, i);
        return f;
    }
};

class ForcePython : public Force {
private:
    py::object pyforce;
    
public:
    ForcePython(py::object _pyforce) : pyforce(_pyforce) {}
    
    void pre_compute(System* system) {
        ExaDisNet disnet(system);
        pyforce.attr("PreCompute")(disnet);
    }
    
    void compute(System* system, bool zero=true) {
        ExaDisNet disnet(system);
        // This can only be called from within ExaDiS modules 
        // (integration, topology, etc.) so we assume that pre_compute()
        // has been called before and there is no need to call it again.
        pyforce.attr("NodeForce")(disnet, false);
    }
    
    Vec3 node_force(System* system, const int& i) {
        SerialDisNet* net = system->get_serial_network();
        NodeTag tag = net->nodes[i].tag;
        ExaDisNet disnet(system);
        return pyforce.attr("OneNodeForce")(disnet, tag).cast<Vec3>();
    }
    
    ~ForcePython() {}
    const char* name() { return "ForcePython"; }
};


/*---------------------------------------------------------------------------
 *
 *    Mobility binding
 *
 *-------------------------------------------------------------------------*/
struct MobilityBind {
    Mobility* mobility = nullptr;
    Params params;
    MobilityBind() {}
    MobilityBind(Mobility* _mobility, Params _params) : 
    mobility(_mobility), params(_params) {}
    std::string name() { return std::string(mobility->name()); }
    void compute(SystemBind& sysbind) { mobility->compute(sysbind.system); }
};

class MobilityPython : public Mobility {
private:
    py::object pymobility;
    
public:
    MobilityPython(py::object _pymobility) : pymobility(_pymobility) {
        non_linear = pymobility.attr("non_linear").cast<bool>(); 
    }
    
    void compute(System* system) {
        ExaDisNet disnet(system);
        pymobility.attr("Mobility")(disnet);
    }
    
    Vec3 node_velocity(System *system, const int& i, const Vec3& fi) {
        SerialDisNet* net = system->get_serial_network();
        NodeTag tag = net->nodes[i].tag;
        ExaDisNet disnet(system);
        return pymobility.attr("OneNodeMobility")(disnet, tag, fi).cast<Vec3>();
    }
    
    ~MobilityPython() {}
    const char* name() { return "MobilityPython"; }
};

template<class M>
MobilityBind make_mobility(Params& params, py::dict mobparams)
{
    params.check_params();
    System* system = make_system(new SerialDisNet(), Crystal(params.crystal), params);
    
    Mobility* mobility = new M(system, mobparams);
    
    exadis_delete(system);
    
    return MobilityBind(mobility, params);
}

} } // namespace ExaDiS::pybind

#include "exadis_pybind_registry.h"
#include "mobility_global_list.h"


namespace ExaDiS { namespace pybind {

/*---------------------------------------------------------------------------
 *
 *    Integrator binding
 *
 *-------------------------------------------------------------------------*/
struct IntegratorBind {
    Integrator* integrator = nullptr;
    Params params;
    IntegratorBind() {}
    IntegratorBind(Integrator* _integrator, Params _params) : 
    integrator(_integrator), params(_params) {}
    double integrate(SystemBind& sysbind) { 
        integrator->integrate(sysbind.system);
        // We also compute plastic strain here
        sysbind.system->plastic_strain();
        // We also reset/update glide planes here
        sysbind.system->reset_glide_planes();
        return sysbind.system->realdt;
    }
};


/*---------------------------------------------------------------------------
 *
 *    Collision binding
 *
 *-------------------------------------------------------------------------*/
struct CollisionBind {
    Collision* collision = nullptr;
    Params params;
    CollisionBind() {}
    CollisionBind(Collision* _collision, Params _params) : 
    collision(_collision), params(_params) {}
    void handle(SystemBind& sysbind) { collision->handle(sysbind.system); }
};


/*---------------------------------------------------------------------------
 *
 *    Topology binding
 *
 *-------------------------------------------------------------------------*/
struct TopologyBind {
    Topology* topology = nullptr;
    Params params;
    double neighbor_cutoff = 0.0;
    TopologyBind() {}
    TopologyBind(Topology* _topology, Params _params, double cutoff) : 
    topology(_topology), params(_params), neighbor_cutoff(cutoff) {}
    void handle(SystemBind& sysbind) { topology->handle(sysbind.system); }
};

typedef typename Topology::Params TParams;

template<class F>
Topology* make_topology_parallel(System* system, Force* force, Mobility* mobility, TParams& topolparams)
{
    Topology* topology = nullptr;
    #define X(NAME, ALIAS) \
        if (strcmp(mobility->name(), #NAME) == 0) { \
            topology = new TopologyParallel<F, MobilityType::ALIAS>(system, force, mobility, topolparams); \
        } else
        EXADIS_MOBILITY_GLOBAL_LIST
        {
            ExaDiS_fatal("Error: invalid mobility type = %s for TopologyParallel binding\n", mobility->name());
        }
    #undef X
    return topology;
}

#include "topology_parallel_types.h"


/*---------------------------------------------------------------------------
 *
 *    Remesh binding
 *
 *-------------------------------------------------------------------------*/
struct RemeshBind {
    Remesh* remesh_class = nullptr;
    Params params;
    RemeshBind() {}
    RemeshBind(Remesh* _remesh, Params _params) : 
    remesh_class(_remesh), params(_params) {}
    void remesh(SystemBind& sysbind) { remesh_class->remesh(sysbind.system); }
};


/*---------------------------------------------------------------------------
 *
 *    Cross-slip binding
 *
 *-------------------------------------------------------------------------*/
struct CrossSlipBind {
    CrossSlip* crossslip = nullptr;
    Params params;
    double neighbor_cutoff = 0.0;
    CrossSlipBind() {}
    CrossSlipBind(CrossSlip* _crossslip, Params _params, double cutoff) : 
    crossslip(_crossslip), params(_params), neighbor_cutoff(cutoff) {}
    void handle(SystemBind& sysbind) { crossslip->handle(sysbind.system); }
};


/*---------------------------------------------------------------------------
 *
 *    Driver binding
 *
 *-------------------------------------------------------------------------*/
class Driver : public ExaDiSApp {
public:
    Driver() : ExaDiSApp() {}
    Driver(const SystemBind& sysbind) : ExaDiSApp() {
        set_system_driver(sysbind);
    }
    virtual void set_system_driver(const SystemBind& sysbind) {
        system = sysbind.system;
    }
    virtual void set_modules_driver(ForceBind& forcebind, MobilityBind& mobbind,
                            IntegratorBind& integratorbind, CollisionBind& collisionbind,
                            TopologyBind& topolbind, RemeshBind& remeshbind,
                            CrossSlipBind& crossslipbind) {
        force = forcebind.force;
        mobility = mobbind.mobility;
        integrator = integratorbind.integrator;
        collision = collisionbind.collision;
        topology = topolbind.topology;
        remesh = remeshbind.remesh_class;
        crossslip = crossslipbind.crossslip;
    }
    virtual py::dict update_state(py::dict& state) {
        double* ptr;
        // edir
        py::array_t<double> pyedir(3);
        ptr = static_cast<double*>(pyedir.request().ptr);
        ptr[0] = edir.x; ptr[1] = edir.y; ptr[2] = edir.z;
        state["edir"] = pyedir;
        // applied_stress
        py::array_t<double> pystress(6);
        ptr = static_cast<double*>(pystress.request().ptr);
        ptr[0] = system->extstress.xx(); ptr[1] = system->extstress.yy(); ptr[2] = system->extstress.zz();
        ptr[3] = system->extstress.yz(); ptr[4] = system->extstress.xz(); ptr[5] = system->extstress.xy();
        state["applied_stress"] = pystress;
        // stress / strain / density
        state["strain"] = strain;
        state["stress"] = stress;
        state["density"] = system->density;
        py::array_t<double> pyEtot(6);
        ptr = static_cast<double*>(pyEtot.request().ptr);
        ptr[0] = Etot.xx(); ptr[1] = Etot.yy(); ptr[2] = Etot.zz();
        ptr[3] = Etot.yz(); ptr[4] = Etot.xz(); ptr[5] = Etot.xy();
        state["Etot"] = pyEtot;
        // time
        state["dt"] = system->realdt;
        state["time"] = tottime;
        state["istep"] = istep;
        return state;
    }
    virtual py::dict read_restart_driver(py::dict& state, std::string restartfile) {
        // read restart
        ExaDiSApp::read_restart(restartfile);
        // set crystal orientation
        py::buffer_info buffer(
            &system->crystal.R, sizeof(double), py::format_descriptor<double>::format(),
            2, {3, 3}, {sizeof(double) * 3, sizeof(double)}
        );
        state["Rorient"] = py::array(buffer);
        // update dictionary
        update_state(state);
        // replace with dummy system so that we don't delete 
        // the original system object upon destruction
        system = make_system(new SerialDisNet());
        return state;
    }
};


} } // namespace ExaDiS::pybind

#endif
