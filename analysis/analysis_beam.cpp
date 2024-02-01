#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include <iostream>

double vector_mean(const std::vector<Double_t>& myVector) {
    // Calculate mean
    double sum = 0.0;
    for (const auto& value : myVector) {
        sum += value;
    }
    double mean = sum / myVector.size();
    return mean;
}

double vector_rms(const std::vector<Double_t>& myVector) {
    double mean = vector_mean(myVector);
    // Calculate root mean square (RMS)
    double sumSquares = 0.0;
    for (const auto& value : myVector) {
        sumSquares += pow(value - mean, 2);
    }
    double rms = sqrt(sumSquares / myVector.size());
    return rms;
}

int analyse_gun() {
    // Open the ROOT file
    TFile* file = TFile::Open("build/data_with_2x50.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening ROOT file" << std::endl;
        return 1;
    }

    // Access the "gun" TTree
    TTree* gunTree = dynamic_cast<TTree*>(file->Get("gun"));
    if (!gunTree) {
        std::cerr << "Error accessing 'gun' TTree" << std::endl;
        file->Close();
        return 1;
    }

    // Variables to store tree data
    Double_t p, x, y, xp, yp, px, py, pz;

    // Set branch addresses
    gunTree->SetBranchAddress("fMom", &p);
    gunTree->SetBranchAddress("fPosX", &x);
    gunTree->SetBranchAddress("fPosY", &y);
    gunTree->SetBranchAddress("fAngX", &xp);
    gunTree->SetBranchAddress("fAngY", &yp);
    gunTree->SetBranchAddress("fMomX", &px);
    gunTree->SetBranchAddress("fMomY", &py);
    gunTree->SetBranchAddress("fMomZ", &pz);

    // Create 2D histograms for x-x' and y-y' scatter plots
    TH2F* scatterX = new TH2F("scatterX", "Scatter plot x-x'; x [mm]; x' [mrad]", 2000, -20, 20, 1000, -50, 50);
    TH2F* scatterY = new TH2F("scatterY", "Scatter plot y-y'; y [mm]; y' [mrad]", 2000, -20, 20, 1000, -50, 50);

    std::vector<Double_t> x_v, xp_v, y_v, yp_v;

    // Loop over events
    Long64_t numEvents = gunTree->GetEntries();
    for (Long64_t i = 0; i < numEvents; ++i) {
        gunTree->GetEntry(i);

        // Fill the vectors
        x_v.push_back(x);
        xp_v.push_back(px/p * 1000);
        y_v.push_back(y);
        yp_v.push_back(py/p * 1000);

        // Fill the scatter plots
        scatterX->Fill(x, px/p*1000);
        scatterY->Fill(y, py/p*1000);
    }

    // Create a canvas for plotting
    TCanvas* canvas = new TCanvas("canvas", "Scatter Plots", 1000, 400);
    canvas->Divide(2, 1);

    // Draw x-x' scatter plot
    canvas->cd(1);
    scatterX->Draw("colz");

    // Draw y-y' scatter plot
    canvas->cd(2);
    scatterY->Draw("colz");

    std::cout << scatterX->GetRMS(1) << " " << vector_rms(x_v) << std::endl;
    std::cout << scatterX->GetRMS(2) << " " << vector_rms(xp_v) << std::endl;
    std::cout << scatterY->GetRMS(1) << " " << vector_rms(y_v) << std::endl;
    std::cout << scatterY->GetRMS(2) << " " << vector_rms(yp_v) << std::endl;

    // Calculate rms emittance
    Double_t emittanceX = scatterX->GetRMS(1) * scatterX->GetRMS(2) /1000;
    Double_t emittanceY = scatterY->GetRMS(1) * scatterY->GetRMS(2) /1000;

    std::cout << "RMS Emittance X: " << emittanceX << " mm*mrad" << std::endl;
    std::cout << "RMS Emittance Y: " << emittanceY << " mm*mrad" << std::endl;

    // Calculate rms emittance
    emittanceX = vector_rms(x_v)*vector_rms(xp_v) /1000;
    emittanceY = vector_rms(y_v)*vector_rms(yp_v) /1000;

    std::cout << "RMS Emittance X: " << emittanceX << " mm*mrad" << std::endl;
    std::cout << "RMS Emittance Y: " << emittanceY << " mm*mrad" << std::endl;

    new TCanvas;
    gunTree->Draw("fMom");

    new TCanvas;
    gunTree->Draw("fPosX:fPosY");

    return 0;
}


int analyse_VD(int req_VDNo) {
    // size of the scint: fPosX > -10 && fPosX < 10 && fPosY > -10 && fPosY < 10
    TCut cutCondition = Form("fVDNo == %d && fInOut == -1 && fParticleID == -13 && fPosX > -10 && fPosX < 10 && fPosY > -10 && fPosY < 10", req_VDNo);

    // Open the ROOT file
    TFile* file = TFile::Open("build/data_with_2x100.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening ROOT file" << std::endl;
        return 1;
    }

    // Access the "VD" TTree
    TTree* VDTree = dynamic_cast<TTree*>(file->Get("VD"));
    if (!VDTree) {
        std::cerr << "Error accessing 'VD' TTree" << std::endl;
        file->Close();
        return 1;
    }

    // Variables to store tree data
    Double_t p, x, y, xp, yp, px, py, pz;
    Int_t VDNo, InOut, ParticleID;

    // Set branch addresses
    VDTree->SetBranchAddress("fParticleID", &ParticleID);
    VDTree->SetBranchAddress("fMom", &p);
    VDTree->SetBranchAddress("fPosX", &x);
    VDTree->SetBranchAddress("fPosY", &y);
    VDTree->SetBranchAddress("fMomX", &px);
    VDTree->SetBranchAddress("fMomY", &py);
    VDTree->SetBranchAddress("fMomZ", &pz);
    VDTree->SetBranchAddress("fVDNo", &VDNo);
    VDTree->SetBranchAddress("fInOut", &InOut);

    // Create 2D histograms for x-x' and y-y' scatter plots
    TH2F* scatterX = new TH2F("scatterX", "Scatter plot x-x'; x [mm]; x' [mrad]", 2000, -20, 20, 1000, -100, 100);
    TH2F* scatterY = new TH2F("scatterY", "Scatter plot y-y'; y [mm]; y' [mrad]", 2000, -20, 20, 1000, -100, 100);

    std::vector<Double_t> x_v, xp_v, y_v, yp_v;

    // Loop over events
    Long64_t numEvents = VDTree->GetEntries();
    for (Long64_t i = 0; i < numEvents; ++i) {
        VDTree->GetEntry(i);
        if(VDNo == req_VDNo && InOut == -1 && ParticleID == -13 && x > -10 && x < 10 && y > -10 && y < 10){
            // Fill the vectors
            x_v.push_back(x);
            xp_v.push_back(px/p * 1000);
            y_v.push_back(y);
            yp_v.push_back(py/p * 1000);

            // Fill the scatter plots
            scatterX->Fill(x, px/p * 1000);
            scatterY->Fill(y, py/p * 1000);
        }
    }

    // Create a canvas for plotting
    TCanvas* canvas = new TCanvas("canvas", "Scatter Plots", 1000, 400);
    canvas->Divide(2, 1);

    // Draw x-x' scatter plot
    canvas->cd(1);
    scatterX->Draw("colz");

    // Draw y-y' scatter plot
    canvas->cd(2);
    scatterY->Draw("colz");

    std::cout << scatterX->GetRMS(1) << " " << vector_rms(x_v) << std::endl;
    std::cout << scatterX->GetRMS(2) << " " << vector_rms(xp_v) << std::endl;
    std::cout << scatterY->GetRMS(1) << " " << vector_rms(y_v) << std::endl;
    std::cout << scatterY->GetRMS(2) << " " << vector_rms(yp_v) << std::endl;

    // Calculate rms emittance
    Double_t emittanceX = scatterX->GetRMS(1) * scatterX->GetRMS(2) /1000;
    Double_t emittanceY = scatterY->GetRMS(1) * scatterY->GetRMS(2) /1000;

    std::cout << "RMS Emittance X: " << emittanceX << " mm*mrad" << std::endl;
    std::cout << "RMS Emittance Y: " << emittanceY << " mm*mrad" << std::endl;

    // Calculate rms emittance
    emittanceX = vector_rms(x_v)*vector_rms(xp_v) /1000;
    emittanceY = vector_rms(y_v)*vector_rms(yp_v) /1000;

    std::cout << "RMS Emittance X: " << emittanceX << " mm*mrad" << std::endl;
    std::cout << "RMS Emittance Y: " << emittanceY << " mm*mrad" << std::endl;

    // Create a canvas for plotting
    TCanvas* canvas2 = new TCanvas("canvas2", "Mom and T", 1000, 400);
    canvas2->Divide(2, 1);

    // Draw x-x' scatter plot
    canvas2->cd(1);
    VDTree->Draw("fMom",cutCondition,"");

    // Draw y-y' scatter plot
    canvas2->cd(2);
    VDTree->Draw("fVDTime",cutCondition,"");

    // Draw y-y' scatter plot
    new TCanvas;
    VDTree->Draw("fPosX:fPosY",cutCondition,"colz");

    return 0;
}

void TOF(){
    // Open the ROOT file
    TFile* file = TFile::Open("build/data.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening ROOT file" << std::endl;
    }

    // Access the "VD" TTree
    TTree* VDTree = dynamic_cast<TTree*>(file->Get("VD"));
    if (!VDTree) {
        std::cerr << "Error accessing 'VD' TTree" << std::endl;
        file->Close();
    }

    // Variables to store tree data
    Double_t p, x, y, xp, yp, px, py, pz, fVDTime;
    Int_t fEvent, VDNo, InOut, ParticleID;

    // Set branch addresses
    VDTree->SetBranchAddress("fParticleID", &ParticleID);
    VDTree->SetBranchAddress("fEvent", &fEvent);
    VDTree->SetBranchAddress("fMom", &p);
    VDTree->SetBranchAddress("fPosX", &x);
    VDTree->SetBranchAddress("fPosY", &y);
    VDTree->SetBranchAddress("fMomX", &px);
    VDTree->SetBranchAddress("fMomY", &py);
    VDTree->SetBranchAddress("fMomZ", &pz);
    VDTree->SetBranchAddress("fVDTime", &fVDTime);
    VDTree->SetBranchAddress("fVDNo", &VDNo);
    VDTree->SetBranchAddress("fInOut", &InOut);

    TH1F* h_tof = new TH1F("h_tof", "h_tof", 300, 0, 30);
    TH1F* h_momIn = new TH1F("h_momIn", "h_momIn", 300, 0, 30);
    TH1F* h_momOut = new TH1F("h_momOut", "h_momOut", 300, 0, 30);

    // Loop over events
    double t_in, t_out;
    int ev_in, ev_out;
    Long64_t numEvents = VDTree->GetEntries();
    for (Long64_t i = 0; i < numEvents; ++i) {
        VDTree->GetEntry(i);
        if(VDNo == 0 && InOut == -1 && ParticleID == -13){ // && x > -10 && x < 10 && y > -10 && y < 10
            ev_in = fEvent;
            t_in = fVDTime;
            h_momIn->Fill(p);
        } 
        if(VDNo == 2 && InOut == -1 && ParticleID == -13){
            ev_out = fEvent;
            t_out = fVDTime;            
            h_momOut->Fill(p);
            if(ev_in == ev_out) h_tof->Fill(t_out-t_in);
        } 
    }
    new TCanvas;
    h_tof->Draw();
    new TCanvas;
    h_momIn->Draw();
    new TCanvas;
    h_momOut->Draw();
}