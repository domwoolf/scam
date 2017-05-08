// SCAM (Soil, Carbon and Microbes) Model ver 1.00
// Author Dominic Woolf
// d.woolf@cornell.edu
//
// Change log:
// Made carbon flux to stable 1st order with mineral surface area (independent of MB)
// Microbial growth rate 1st order with microbial biomass and DOC, limited by competition with mineral adsorption.
// Partitioning between CO2 and microbial growth by CUE (carbon use efficiency)
// CUE linear with temperature
// mineral sorption rate constant linear with clay

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Calculate acquired topsoil moisture deficit, based on pet (potential evapotranspiration) and precip (precipitation)
// Only used when SCAM is called with use_atsmd = TRUE
// otherwise, SCAM uses h2O and sat vectors (soil water content and saturation capacity)
double calc_atsmd (double prec, double pet, double max_tsmd, int cover, double init_tsmd){
    double atsmd = 0.0;
    double excess = prec - pet;
    if (cover) atsmd = std::max(init_tsmd + excess, max_tsmd);
    else if (init_tsmd < max_tsmd/1.8) atsmd = init_tsmd + std::max(excess, 0.0);
    return(std::min(atsmd, 0.0));
}

double calc_clayfact (double clay, double mclay) {
    return(mclay*clay + 1 - 23*mclay);
}

double calc_cue (double temp, double cue_0, double mcue) {
    //http://onlinelibrary.wiley.com.proxy.library.cornell.edu/doi/10.1890/15-2110.1/full
    //http://onlinelibrary.wiley.com/doi/10.1111/ele.12113/full
    //http://link.springer.com/article/10.1007/s10533-013-9948-8
    return(cue_0 - (temp-15)*mcue); // default cue_0 = 0.31, mcue=0.005
    //return (0.6 - std::exp(-0.05 * std::exp(0.11*temp) - 0.7));
}

// [[Rcpp::export]]
List scam(
    DataFrame soc_data, 
    double mic_fact, 
    double mic_offset,
    double kdissolution, 
    double kdepoly, 
    double kdeath_and_exudates, 
    double kdesorb, 
    double ksorb, 
    double kmicrobial_uptake,
    double cue_0,
    double mcue,
    double mclay,
    double clay,
    double max_tsmd = 0.0, 
    bool use_atsmd = 0,
    bool monthly = 0) {

    // extract input dataframe columns into C++ vectors
    IntegerVector day = soc_data["day"];
    IntegerVector month = soc_data["month"];
    NumericVector dpm = soc_data["dpm"];
    NumericVector rpm = soc_data["rpm"];
    NumericVector doc = soc_data["doc"];
    NumericVector bio = soc_data["bio"];
    NumericVector sta = soc_data["sta"];
    NumericVector soc = soc_data["soc"];
    NumericVector add_dpm = soc_data["added.dpm"];
    NumericVector add_rpm = soc_data["added.rpm"];
    NumericVector add_doc = soc_data["added.doc"];
    NumericVector add_bio = soc_data["added.bio"];
    NumericVector add_sta = soc_data["added.hum"];
    IntegerVector cover = soc_data["cover"];
    NumericVector atsmd = soc_data["atsmd"];
    NumericVector temp = soc_data["temp"];
    NumericVector precip = soc_data["precip"];
    NumericVector pet = soc_data["pet"];
    NumericVector sat = soc_data["sat"];
    NumericVector h2o = soc_data["h2o"];
    NumericVector a = soc_data["a"];
    NumericVector c = soc_data["c"];
    NumericVector add_d13c = soc_data["added.d13c"];
    NumericVector soc_d13c = soc_data["soc.d13c"];
    
    // define new vectors
    int n = day.size();
    NumericVector dec_dpm(n), dec_rpm(n), dec_doc(n), dec_bio(n), dec_sta(n), dec_tot(n);
    NumericVector sorption(n), microbial_uptake(n), growth(n);
    NumericVector min_doc(n), min_cum(n);
    NumericVector b(n), mic(n), cue(n);
    NumericVector dpm_d13c(n), rpm_d13c(n), doc_d13c(n), bio_d13c(n), sta_d13c(n), co2_d13c(n);    
    
    double clayfact = calc_clayfact(clay, mclay);
    double ksorb_altered, kmicrobial_uptake_altered, fsorb, fmic, kdoc;

    for (unsigned int i = 0; i < n; i++) { 
        if (i > 0) {  // first month has pool sizes and delta 13C values already initialized. Calculate only for subsequent months
            // pools
            dpm[i] = dpm[i-1] - dec_dpm[i-1] + add_dpm[i-1];
            rpm[i] = rpm[i-1] - dec_rpm[i-1] + add_rpm[i-1];
            doc[i] = doc[i-1] - dec_doc[i-1] + add_doc[i-1] + dec_dpm[i-1] + dec_rpm[i-1] + dec_bio[i-1] + dec_sta[i-1];
            bio[i] = bio[i-1] - dec_bio[i-1] + add_bio[i-1] + growth[i-1];
            sta[i] = sta[i-1] - dec_sta[i-1] + add_sta[i-1] + sorption[i-1];
            soc[i] = dpm[i] + rpm[i] + doc[i] + bio[i] + sta[i];
            
            //d13C of soil carbon pools    
            dpm_d13c[i] = ((dpm[i-1] - dec_dpm[i-1])*dpm_d13c[i-1] + add_dpm[i-1]*add_d13c[i-1]) / dpm[i];
            rpm_d13c[i] = ((rpm[i-1] - dec_rpm[i-1])*rpm_d13c[i-1] + add_rpm[i-1]*add_d13c[i-1]) / rpm[i];
            doc_d13c[i] = ((doc[i-1] - dec_doc[i-1])*doc_d13c[i-1] + add_doc[i-1]*add_d13c[i-1] + dec_dpm[i-1] * dpm_d13c[i-1] + dec_rpm[i-1] * rpm_d13c[i-1] + dec_bio[i-1] * bio_d13c[i-1] + dec_sta[i-1] * sta_d13c[i-1]) / doc[i];
            bio_d13c[i] = ((bio[i-1] - dec_bio[i-1]) * bio_d13c[i-1] + growth[i-1] * doc_d13c[i-1]) / bio[i];
            sta_d13c[i] = ((sta[i-1] - dec_sta[i-1]) * sta_d13c[i-1] + add_sta[i-1] * add_d13c[i-1] + sorption[i-1] * doc_d13c[i-1]) / sta[i];
        }
        
        // rate modifiers for moisture (b) and for microbial biomass (mic)
        if (use_atsmd) {
            atsmd[i] = (i == 0) ? 0 : calc_atsmd(precip[i], pet[i], max_tsmd, cover[i], atsmd[i-1]);
            b[i] = (atsmd[i] > 0.444*max_tsmd) ? 1 : (0.2 + 0.8*(max_tsmd-atsmd[i]) / (0.556*max_tsmd));
        }
        else b[i] =  (h2o[i] > 0.556*sat[i]) ? 1 : 0.2 + 0.8 * h2o[i] / (0.556*sat[i]);
        mic[i] = mic_fact*pow(bio[i], mic_offset);
        
        // partitioning of doc between competing processes of sorption and microbial uptake
        ksorb_altered = ksorb * clayfact * a[i] * b[i] * c[i]; // sorption rate multiplied by rate modifiers
        kmicrobial_uptake_altered = kmicrobial_uptake * mic[i] * a[i] * b[i] * c[i]; // microbial uptake multiplied by modifiers
        fsorb = ksorb_altered / (ksorb_altered + kmicrobial_uptake_altered); // fraction of doc turnover sorbed
        fmic = 1 - fsorb; // fraction of doc turnover taken up by microbes
        kdoc = (fsorb * ksorb_altered) + (fmic * kmicrobial_uptake_altered); // rate constant for doc loss is weighted mean of sorption and microbial uptake
        dec_doc[i] = doc[i] * (1 - std::exp(-kdoc));
        sorption[i] = dec_doc[i] * fsorb;
        microbial_uptake[i] = dec_doc[i] * fmic;
        
        // partitioning of microbial uptake between growth and respiration
        cue[i] = calc_cue(temp[i], cue_0, mcue);
        growth[i] = cue[i] * microbial_uptake[i];
        min_doc[i] = (1- cue[i]) * microbial_uptake[i];         // Mineralization in time period
        min_cum[i] = i==0 ? min_doc[i] : min_cum[i-1] + min_doc[i];   // cumulative mineralization 
        
        // Decomposition of non-doc pools
        dec_dpm[i] = dpm[i] * (1 - std::exp(-kdissolution * a[i] * b[i] * c[i] * mic[i]));
        dec_rpm[i] = rpm[i] * (1 - std::exp(-kdepoly * a[i] * b[i] * c[i] * mic[i]));
        dec_bio[i] = bio[i] * (1 - std::exp(-kdeath_and_exudates * a[i] * b[i] * c[i] * mic[i]));
        dec_sta[i] = sta[i] * (1 - std::exp(-kdesorb * a[i] * b[i] * c[i] * mic[i]));
        dec_tot[i] = dec_dpm[i] + dec_rpm[i] + dec_doc[i] + dec_bio[i] + dec_sta[i];        

        //d13C of respired CO2 and total soil organic carbon 
        co2_d13c[i] = doc_d13c[i];
        soc_d13c[i] = (dpm[i]*dpm_d13c[i] + rpm[i]*rpm_d13c[i] + doc[i]*doc_d13c[i] + bio[i]*bio_d13c[i] + sta[i]*sta_d13c[i]) / soc[i];
    }

    //export values back to R as a list
    List output;
    output["day"] = day;
    output["month"] = month;
    output["dpm"] = dpm;
    output["rpm"] = rpm;
    output["doc"] = doc;
    output["bio"] = bio;
    output["sta"] = sta;
    output["soc"] = soc;
    output["atsmd"] = atsmd;
    output["cue"] = cue;
    output["b"] = b;
    output["mic"] = mic;
    output["sorption"] = sorption;
    output["microbial_uptake"] = microbial_uptake;
    output["growth"] = growth;
    output["decomp.dpm"] = dec_dpm;
    output["decomp.rpm"] = dec_rpm;
    output["decomp.doc"] = dec_doc;
    output["decomp.bio"] = dec_bio;
    output["decomp.sta"] = dec_sta;
    output["decomp.tot"] = dec_tot;
    output["mineralise.doc"] = min_doc;
    output["mineralise.cum"] = min_cum;
    output["dpm.d13c"] = dpm_d13c;
    output["rpm.d13c"] = rpm_d13c;
    output["doc.d13c"] = doc_d13c;
    output["bio.d13c"] = bio_d13c;
    output["sta.d13c"] = sta_d13c;
    output["soc.d13c"] = soc_d13c;
    output["co2.d13c"] = co2_d13c;
    return output;  
}
