[![R-CMD-check-bioc](../../actions/workflows/check-bioc.yml/badge.svg)](../../actions?query=workflow%3AR-CMD-check-bioc)

# Intoduction 

Because of the numerous benefits breastfeeding provides, it’s recommended that infants be exclusively breastfed for the first 6 months with continued breastfeeding for at least 1 year. However lactation is not achieved in all mothers and the physiology of this process remains poorly understood. This is because of questions related to the ethics of performing tests or extracting tissues from lactating mother.

So, an interesting alternative has been proposed – human milk secreted during lactation is a rich source of mammary epithelial cell RNA. About 3–8% of human milk fat globules contain mammary epithelial cell cytoplasmic remnants, including RNA, captured during milk fat globule formation and secretion

Previous to this research, there has been attempt to undestand the human milk fat layer transcriptome with microarray techniques (Maningat et al. 2009). This paper analyses the transcriptome with RNA-seq which, compared to microarray, has highly sensitive detail, can accurately quantify both very low and very high abundance transcripts, as well as detect novel transcripts.

The authors of this paper (Lemay et al. 2013), furthermore, suggest that the constitution of the human milk is not dependent on the days postpartum, which is the measure typically used by lactation researchers to categorize lactation stage. Rather, they state that compared to other factors such as subject, washing protocol, or time of collection, biochemically defined lactation stage (i.e., Na:K) has the largest effect on the transcriptome. So they categorise the samples as: “colostral” if Na:K >= 2.0; and as “transitional” if Na:K < 2.0; “mature” if Na:K < 0.6 .

The three main challanges of this paper are:

*    because of the sample extraction, it is unknown if the milk fat layer yields RNA purely of mammary epithelial cell origin, or if RNA from non-mammary epithelial sources is also present in the milk fat layer obtained from colostrum, transitional, or mature human milk,

*    because of the body heat and the environment of the milk, the RNA might end up being degraded (low RIN values),

*    because this is a study of patients who would agree to participate after just having gone through giving birth and have to fulfill some criteria, the sample size is very small which is something that really influences the information that can be extracted and extrapolated from this sample group - the result of bias of such a small sample group is really apparent in the MA-plots


# Code for analysis

The code is contained it the directory 'vignettes' (file: *IEOprojectAnalysis.Rmd*).

# Built html page

The built html file from this code is located in the directory 'doc' (file: *IEOprojectAnalysis.html*).

# Reference

The data and idea for the approach are taken from the paper:

Lemay, D. G., O. A. Ballard, M. A. Hughes, A. L. Morrow, N. D. Lemay DG, Ballard OA, Hughes MA, Morrow AL, Horseman ND, Nommsen-Rivers LA. RNA sequencing of the human milk fat layer transcriptome reveals distinct gene expression profiles at three stages of lactation. PLoS One. 2013 Jul 5;8(7):e67531. doi: 10.1371/journal.pone.0067531. PMID: 23861770; PMCID: PMC3702532.

