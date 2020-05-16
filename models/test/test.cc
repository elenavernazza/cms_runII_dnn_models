// C++
#include <set>
#include <stdexcept>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <cmath>

// ROOT
#include <Math/VectorUtil.h>
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/PxPyPzM4D.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

// Plugins
#include "cms_hh_tf_inference/inference/interface/inf_wrapper.hh"
#include "cms_hh_proc_interface/processing/interface/feat_comp.hh"
#include "cms_hh_proc_interface/processing/interface/evt_proc.hh"

const double E_MASS  = 0.0005109989; //GeV
const double MU_MASS = 0.1056583715; //GeV

std::string model_dir = "../../src/cms_runII_dnn_models/models/nonres_gluglu/";
std::string data_dir = "/eos/home-k/kandroso/cms-it-hh-bbtautau/anaTuples/";

using LorentzVectorPEP = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;
using LorentzVector    = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;


std::vector<std::string>get_evt_names(const std::map<unsigned long, std::string>& id2name, const std::vector<unsigned long>& ids) {
    /* Match data IDs to aux names */
    
    std::vector<std::string> names(ids.size());
    for (unsigned int i = 0; i < ids.size(); i++) names[i] = id2name.at(ids[i]);
    return names;
}


int jet_cat_lookup(const std::string& jet_cat) {
    if (jet_cat == "2j")            return 0;
    if (jet_cat == "2j0bR_noVBF")   return 1;
    if (jet_cat == "2j1bR_noVBF")   return 2;
    if (jet_cat == "2j2b+R_noVBF")  return 3;
    if (jet_cat == "2j2Lb+B_noVBF") return 4;
    if (jet_cat == "2j1b+_VBFL")    return 5;
    if (jet_cat == "2j1b+_VBF")     return 6;
    throw std::invalid_argument("Unrecognised jet category: " + jet_cat);
    return -1;
}


int region_lookup(const std::string& region) {
    if (region == "OS_Isolated")      return 0;
    if (region == "OS_AntiIsolated")  return 1;
    if (region == "SS_Isolated")      return 2;
    if (region == "SS_AntiIsolated")  return 3;
    if (region == "SS_LooseIsolated") return 4;
    throw std::invalid_argument("Unrecognised region: " + region);
    return -1;
}

void sample_lookup(const std::string& sample, int& sample_id, Spin& spin, float& klambda, float& res_mass) {
    spin = nonres;
    res_mass = 125;
    klambda = 1;
    
    if (sample.find("GluGluSignal") != std::string::npos) { 
        if (sample.find("NonRes") != std::string::npos) {
            sample_id = -12;
            try {
                klambda = std::stof(sample.substr(sample.find("_kl")+3));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_kl")+3) << "\n";
                assert(false);
            }
        } else if (sample.find("Radion") != std::string::npos) {
            spin = radion;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -13;
            } else if (res_mass <= 600) {
                sample_id = -14;
            } else {
                sample_id = -15;
            }
        } else if (sample.find("Graviton") != std::string::npos) {
            spin = graviton;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -16;
            } else if (res_mass <= 600) {
                sample_id = -17;
            } else {
                sample_id = -18;
            }
        }
    } else if (sample.find("VBFSignal") != std::string::npos) {
        if (sample.find("NonRes") != std::string::npos) {
            sample_id = -19;
        } else if (sample.find("Radion") != std::string::npos) {
            spin = radion;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -20;
            } else if (res_mass <= 600) {
                sample_id = -21;
            } else {
                sample_id = -22;
            }
        } else if (sample.find("Graviton") != std::string::npos) {
            spin = graviton;
            try {
                res_mass = std::stof(sample.substr(sample.find("_M")+2));
            } catch (...) {
                std::cout << "Error in sample " << sample << " attempting to parse " << sample.substr(sample.find("_M")+2) << "\n";
                assert(false);
            }
            if (res_mass <= 400) {
                sample_id = -23;
            } else if (res_mass <= 600) {
                sample_id = -24;
            } else {
                sample_id = -25;
            }
        }
    } else if (sample.find("Data") != std::string::npos) {
        sample_id = 0;
    } else if (sample.find("TT") != std::string::npos) {
        sample_id = 1;
    } else if (sample.find("ttH") != std::string::npos) {
        sample_id = 2;
    } else if (sample.find("DY") != std::string::npos) {
        sample_id = 3;
    } else if (sample.find("Wjets") != std::string::npos) {
        sample_id = 4;
    } else if (sample.find("SM_Higgs") != std::string::npos) {
        sample_id = 5;
    } else if (sample.find("VH") != std::string::npos) {
        sample_id = 6;
    } else if (sample.find("VVV") != std::string::npos) {
        sample_id = 7;
    } else if (sample.find("EWK") != std::string::npos) {
        sample_id = 8;
    } else if (sample.find("VV") != std::string::npos) {
        sample_id = 9;
    } else if (sample.find("ST") != std::string::npos) {
        sample_id = 10;
    } else if (sample.find("ttV") != std::string::npos) {
        sample_id = 11;
    } else{
        throw std::invalid_argument("Unrecognised sample: " + sample);
    }
}


int sample2class_lookup(const int& sample) {
    if (sample < 0)  return 1;   // Signal
    if (sample == 0) return -1;  // Collider data
    return 0;                    // Background
}


void extract_flags(const std::vector<std::string>& name, int& sample, int& region, bool& syst_unc, bool& scale, int& jet_cat, int& cut, int& class_id,
                   Spin& spin, float& klambda, float& res_mass, bool& is_boosted) {
    /*
    Extract event flags from name 
    Example: "2j/NoCuts/SS_AntiIsolated/None/Central/DY_MC_M-10-50"
    */

    std::string val;
    int tmp;
    jet_cat = -1;
    cut = -1;
    for (unsigned int n = 0; n < name.size(); n++) {
        std::istringstream iss(name[n]);
        int i = 0;
        while (std::getline(iss, val, '/')) {   
            if (i == 0) {
                tmp = jet_cat_lookup(val);
                if (tmp > jet_cat) jet_cat = tmp;
                is_boosted = (jet_cat == 5);
            } else if (i == 1) {
                tmp = cut_lookup(val);
                if (tmp > cut) cut = tmp;
            } else if (i == 2 && n == 0) {
                region = region_lookup(val);
            } else if (i == 3 && n == 0) {
                syst_unc = (val == "None");
            } else if (i == 4 && n == 0) {
                scale = (val == "Central");
            } else if (i == 5 && n == 0) {
                sample_lookup(val, sample, spin, klambda, res_mass);
            }
            i++;
        }
        class_id = sample2class_lookup(sample);
    }
}


std::map<unsigned long, std::string> build_id_map(TFile* in_file) {
    TTreeReader aux_reader("aux", in_file);
    TTreeReaderValue<std::vector<unsigned long>> rv_aux_id(aux_reader, "dataIds");
    TTreeReaderValue<std::vector<std::string>> rv_aux_name(aux_reader, "dataId_names");
    std::vector<unsigned long> ids;
    std::vector<std::string> names;
    while (aux_reader.Next()) {
        ids   = *rv_aux_id;
        names = *rv_aux_name;
    }
    std::map<unsigned long, std::string> id2name;
    for (unsigned int i = 0; i < ids.size(); i++) id2name[ids[i]] = names[i];
    return id2name;
}


std::vector<std::string> get_requested(std::string feat_file) {
    std::ifstream infile(feat_file);
    std::cout << "Reading features from file: " << feat_file << "\nFeatures:";
    std::string line;
    std::vector<std::string> requested;
    while (std::getline(infile, line)) {
        std::cout << line << " ";
        requested.push_back(line);
    }
    std::cout << "\n";
    infile.close();
    return requested;
}


bool run_test_loop(std::string fname, InfWrapper wrapper, const int& n, const std::string& feat_file) {
    std::cout << "Reading data from file: " << fname << "\n";
    TFile* in_file = TFile::Open(fname.c_str());
    TTreeReader reader("muTau", in_file);
    
    std::vector<std::string> requested = get_requested(feat_file);

    EvtProc evt_proc(false, requested, true);

    // Enums
    Channel e_channel(muTau);
    std::string channel = "muTau";
    Year e_year(y18);
    Spin spin(nonres);
    float klambda;
    float res_mass;

    // Meta info
    std::cout << "Extracting auxiliary data...";
    TTreeReaderValue<unsigned long long> rv_evt(reader, "evt");
    TTreeReaderValue<std::vector<unsigned long>> rv_id(reader, "dataIds");
    std::map<unsigned long, std::string> id2name = build_id_map(in_file);
    std::cout << " Extracted\n";
    std::vector<std::string> names;
    int sample, region, jet_cat, cut, n_vbf, class_id;
    unsigned long long int evt;
    bool scale, syst_unc, svfit_conv, hh_kinfit_conv;
    std::vector<unsigned long> ids;

    // HL feats
    TTreeReaderValue<float> rv_kinfit_mass(reader, "kinFit_m");
    TTreeReaderValue<float> rv_kinfit_chi2(reader, "kinFit_chi2");
    TTreeReaderValue<float> rv_mt2(reader, "MT2");
    // TTreeReaderValue<float> rv_mt_tot(reader, "mt_tot");
    // TTreeReaderValue<float> rv_top_1_mass(reader, "mass_top1");
    // TTreeReaderValue<float> rv_top_2_mass(reader, "mass_top2");
    // TTreeReaderValue<float> rv_p_zetavisible(reader, "p_zetavisible");
    // TTreeReaderValue<float> rv_p_zeta(reader, "p_zeta");
    float kinfit_mass, kinfit_chi2, mt2, mt_tot, top_1_mass, top_2_mass, p_zetavisible, p_zeta;

    // Tagging
    TTreeReaderValue<float> rv_b_1_csv(reader, "b1_DeepFlavour");
    TTreeReaderValue<float> rv_b_2_csv(reader, "b2_DeepFlavour");
    float b_1_csv, b_2_csv;
    bool is_boosted;

    // SVFit feats
    TTreeReaderValue<float> rv_svfit_pT(reader, "SVfit_pt");
    TTreeReaderValue<float> rv_svfit_eta(reader, "SVfit_eta");
    TTreeReaderValue<float> rv_svfit_phi(reader, "SVfit_phi");
    TTreeReaderValue<float> rv_svfit_mass(reader, "SVfit_m");
    LorentzVectorPEP pep_svfit;
    LorentzVector svfit;

    // l1 feats
    TTreeReaderValue<float> rv_l_1_pT(reader, "tau1_pt");
    TTreeReaderValue<float> rv_l_1_eta(reader, "tau1_eta");
    TTreeReaderValue<float> rv_l_1_phi(reader, "tau1_phi");
    TTreeReaderValue<float> rv_l_1_mass(reader, "tau1_m");
    // TTreeReaderValue<float> rv_l_1_mt(reader, "mt_1");
    float l_1_mass, l_1_mt;
    LorentzVectorPEP pep_l_1;
    LorentzVector l_1;


    // l2 feats
    TTreeReaderValue<float> rv_l_2_pT(reader, "tau2_pt");
    TTreeReaderValue<float> rv_l_2_eta(reader, "tau2_eta");
    TTreeReaderValue<float> rv_l_2_phi(reader, "tau2_phi");
    TTreeReaderValue<float> rv_l_2_mass(reader, "tau2_m");
    // TTreeReaderValue<float> rv_l_2_mt(reader, "mt_2");
    float l_2_mt;
    LorentzVectorPEP pep_l_2;\
    LorentzVector l_2;

    // MET feats
    TTreeReaderValue<float> rv_met_pT(reader, "MET_pt");
    TTreeReaderValue<float> rv_met_phi(reader, "MET_phi");
    LorentzVectorPEP pep_met;
    LorentzVector met;

    // b1 feats
    TTreeReaderValue<float> rv_b_1_pT(reader, "b1_pt");
    TTreeReaderValue<float> rv_b_1_eta(reader, "b1_eta");
    TTreeReaderValue<float> rv_b_1_phi(reader, "b1_phi");
    TTreeReaderValue<float> rv_b_1_mass(reader, "b1_m");
    TTreeReaderValue<float> rv_b_1_hhbtag(reader, "b1_HHbtag");
    TTreeReaderValue<float> rv_b_1_cvsl(reader, "b1_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_b_1_cvsb(reader, "b1_DeepFlavour_CvsB");
    float b_1_hhbtag, b_1_cvsl, b_1_cvsb;
    LorentzVectorPEP pep_b_1;
    LorentzVector b_1;

    // b2 feats
    TTreeReaderValue<float> rv_b_2_pT(reader, "b2_pt");
    TTreeReaderValue<float> rv_b_2_eta(reader, "b2_eta");
    TTreeReaderValue<float> rv_b_2_phi(reader, "b2_phi");
    TTreeReaderValue<float> rv_b_2_mass(reader, "b2_m");
    TTreeReaderValue<float> rv_b_2_hhbtag(reader, "b2_HHbtag");
    TTreeReaderValue<float> rv_b_2_cvsl(reader, "b2_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_b_2_cvsb(reader, "b2_DeepFlavour_CvsB");
    float b_2_hhbtag, b_2_cvsl, b_2_cvsb;
    LorentzVectorPEP pep_b_2;
    LorentzVector b_2;

    // vbf1 feats
    TTreeReaderValue<float> rv_vbf_1_pT(reader, "VBF1_pt");
    TTreeReaderValue<float> rv_vbf_1_eta(reader, "VBF1_eta");
    TTreeReaderValue<float> rv_vbf_1_phi(reader, "VBF1_phi");
    TTreeReaderValue<float> rv_vbf_1_mass(reader, "VBF1_m");
    TTreeReaderValue<float> rv_vbf_1_hhbtag(reader, "VBF1_HHbtag");
    TTreeReaderValue<float> rv_vbf_1_cvsl(reader, "VBF1_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_vbf_1_cvsb(reader, "VBF1_DeepFlavour_CvsB");
    float vbf_1_hhbtag, vbf_1_cvsl, vbf_1_cvsb;
    LorentzVectorPEP pep_vbf_1;
    LorentzVector vbf_1;

    // vbf2 feats
    TTreeReaderValue<float> rv_vbf_2_pT(reader, "VBF2_pt");
    TTreeReaderValue<float> rv_vbf_2_eta(reader, "VBF2_eta");
    TTreeReaderValue<float> rv_vbf_2_phi(reader, "VBF2_phi");
    TTreeReaderValue<float> rv_vbf_2_mass(reader, "VBF2_m");
    TTreeReaderValue<float> rv_vbf_2_hhbtag(reader, "VBF2_HHbtag");
    TTreeReaderValue<float> rv_vbf_2_cvsl(reader, "VBF2_DeepFlavour_CvsL");
    TTreeReaderValue<float> rv_vbf_2_cvsb(reader, "VBF2_DeepFlavour_CvsB");
    float vbf_2_hhbtag, vbf_2_cvsl, vbf_2_cvsb;
    LorentzVectorPEP pep_vbf_2;
    LorentzVector vbf_2;

    std::vector<float> feat_vals;
    float pred;
    std::cout << "\tprepared.\nBeginning loop.\n";

    long int c_event(0), n_tot_events(reader.GetEntries(true));
    while (reader.Next()) {
        c_event++;
        if (c_event%1000 == 0) std::cout << c_event << " / " << n_tot_events << "\n";
        ids = *rv_id;

        names = get_evt_names(id2name, ids);
        extract_flags(names, sample, region, syst_unc, scale, jet_cat, cut, class_id, spin, klambda, res_mass, is_boosted);

        // Load meta
        evt    =  *rv_evt;

        // Load HL feats
        kinfit_mass   = *rv_kinfit_mass;
        kinfit_chi2   = *rv_kinfit_chi2;
        mt2           = *rv_mt2;
        mt_tot        = 0;  // *rv_mt_tot;
        top_1_mass    = 0;  // *rv_top_1_mass;
        top_2_mass    = 0;  // *rv_top_2_mass;
        p_zetavisible = 0;  // *rv_p_zetavisible;
        p_zeta        = 0;  // *rv_p_zeta;
        l_1_mt        = 0;  // *rv_l_1_mt;
        l_2_mt        = 0;  // *rv_l_2_mt;
        b_1_hhbtag    = *rv_b_1_hhbtag;
        b_2_hhbtag    = *rv_b_2_hhbtag;
        vbf_1_hhbtag  = *rv_vbf_1_hhbtag;
        vbf_2_hhbtag  = *rv_vbf_2_hhbtag;
        b_1_cvsl      = *rv_b_1_cvsl;
        b_2_cvsl      = *rv_b_2_cvsl;
        vbf_1_cvsl    = *rv_vbf_1_cvsl;
        vbf_2_cvsl    = *rv_vbf_2_cvsl;
        b_1_cvsb      = *rv_b_1_cvsb;
        b_2_cvsb      = *rv_b_2_cvsb;
        vbf_1_cvsb    = *rv_vbf_1_cvsb;
        vbf_2_cvsb    = *rv_vbf_2_cvsb;

        // Load tagging
        b_1_csv     = *rv_b_1_csv;
        b_2_csv     = *rv_b_2_csv;

        // Load vectors
        pep_svfit.SetCoordinates(*rv_svfit_pT, *rv_svfit_eta, *rv_svfit_phi, *rv_svfit_mass);
        if (channel == "muTau") {  // Fix mass for light leptons
            l_1_mass = MU_MASS;
        } else if (channel == "eTau") {
            l_1_mass = E_MASS;
        } else {
            l_1_mass = *rv_l_1_mass;
        }
        pep_l_1.SetCoordinates(*rv_l_1_pT, *rv_l_1_eta, *rv_l_1_phi, l_1_mass);
        pep_l_2.SetCoordinates(*rv_l_2_pT, *rv_l_2_eta, *rv_l_2_phi, *rv_l_2_mass);
        pep_met.SetCoordinates(*rv_met_pT, 0,           *rv_met_phi, 0);
        pep_b_1.SetCoordinates(*rv_b_1_pT, *rv_b_1_eta, *rv_b_1_phi, *rv_b_1_mass);
        pep_b_2.SetCoordinates(*rv_b_2_pT, *rv_b_2_eta, *rv_b_2_phi, *rv_b_2_mass);
        pep_vbf_1.SetCoordinates(*rv_vbf_1_pT, *rv_vbf_1_eta, *rv_vbf_1_phi, *rv_vbf_1_mass);
        pep_vbf_2.SetCoordinates(*rv_vbf_2_pT, *rv_vbf_2_eta, *rv_vbf_2_phi, *rv_vbf_2_mass);

        svfit.SetCoordinates(pep_svfit.Px(), pep_svfit.Py(), pep_svfit.Pz(), pep_svfit.M());
        l_1.SetCoordinates(pep_l_1.Px(),     pep_l_1.Py(),   pep_l_1.Pz(),   pep_l_1.M());
        l_2.SetCoordinates(pep_l_2.Px(),     pep_l_2.Py(),   pep_l_2.Pz(),   pep_l_2.M());
        met.SetCoordinates(pep_met.Px(),     pep_met.Py(),   0,              0);
        b_1.SetCoordinates(pep_b_1.Px(),     pep_b_1.Py(),   pep_b_1.Pz(),   pep_b_1.M());
        b_2.SetCoordinates(pep_b_2.Px(),     pep_b_2.Py(),   pep_b_2.Pz(),   pep_b_2.M());
        vbf_1.SetCoordinates(pep_vbf_1.Px(), pep_vbf_1.Py(), pep_vbf_1.Pz(), pep_vbf_1.M());
        vbf_2.SetCoordinates(pep_vbf_2.Px(), pep_vbf_2.Py(), pep_vbf_2.Pz(), pep_vbf_2.M());

        // VBF
        n_vbf = 0;
        if (jet_cat == 4) {
            if (*rv_vbf_1_mass != std::numeric_limits<float>::lowest()) n_vbf++;
            if (*rv_vbf_2_mass != std::numeric_limits<float>::lowest()) n_vbf++;
        }

        // Convergence
        svfit_conv     = *rv_svfit_mass > 0;
        hh_kinfit_conv = kinfit_chi2    > 0;

        feat_vals = evt_proc.process_as_vec(b_1, b_2, l_1, l_2, met, svfit, vbf_1, vbf_2, kinfit_mass, kinfit_chi2, mt2, mt_tot, p_zetavisible, p_zeta,
                                            top_1_mass, top_2_mass, l_1_mt, l_2_mt, is_boosted, b_1_csv, b_2_csv, e_channel, e_year,
                                            res_mass, spin, klambda, n_vbf, svfit_conv, hh_kinfit_conv, b_1_hhbtag, b_2_hhbtag, vbf_1_hhbtag, vbf_2_hhbtag,
                                            b_1_cvsl, b_2_cvsl, vbf_1_cvsl, vbf_2_cvsl, b_1_cvsb, b_2_cvsb, vbf_1_cvsb, vbf_2_cvsb, cv, c2v, c3);

        std::cout << "Input features:\n";
        for (unsigned int i=0; i < requested.size(); i++) std::cout << requested[i] << "\t:\t" << feat_vals[i] << "\n";
        pred = wrapper.predict(feat_vals, evt);
        std::cout << "\nEvent " << c_event << " class " << class_id << " prediction " << pred << "\n";

        if (c_event >= n && n > 0) break;
    }

    std::cout << "Loop complete.\n";
    in_file->Close();
    return true;
}

void show_help() {
    /* Show help for input arguments */

    std::cout << "-n : number of events to run, default 1, set negative to run all events in file\n";
    std::cout << "-m : model version to use, default 2020-03-11-0\n";
    std::cout << "-d : data version to use, default 2020-02-14\n";
}

std::map<std::string, std::string> get_options(int argc, char* argv[]) {
    /*Interpret input arguments*/

    std::map<std::string, std::string> options;
    options.insert(std::make_pair("-n", "1")); // number of events
    options.insert(std::make_pair("-m", "2020-03-11-0")); // moddel version
    options.insert(std::make_pair("-d", "2020-02-14")); // moddel version

    if (argc >= 2) { //Check if help was requested
        std::string option(argv[1]);
        if (option == "-h" || option == "--help") {
            show_help();
            options.clear();
            return options;
        }
    }

    for (int i = 1; i < argc; i = i+2) {
        std::string option(argv[i]);
        std::string argument(argv[i+1]);
        if (option == "-h" || option == "--help" || argument == "-h" || argument == "--help") { // Check if help was requested
            show_help();
            options.clear();
            return options;
        }
        options[option] = argument;
    }
    return options;
}

int main(int argc, char *argv[]) {
    std::map<std::string, std::string> options = get_options(argc, argv); // Parse arguments
    if (options.size() == 0) return 1;

    std::cout << "Instantiating wrapper\n";
    InfWrapper wrapper(model_dir+options["-m"]+"/ensemble", 1, true);
    std::cout << "Wrapper instantiated\n";

    std::cout << "\nBeginning test loop for ensemble\n";
    assert(run_test_loop(data_dir+options["-d"]+"/2018_muTau.root", wrapper, std::stoi(options["-n"]), model_dir+options["-m"]+"/features.txt"));
    std::cout << "\nAll tests completed sucessfully\n";
    return 0;
}