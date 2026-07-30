# De-identified-Influenza-Large-Dataset
Influenza A and B surveillance and research data set with sample collection dates starting in 2021 with genomic sequencing information as well as corresponding demographic, geolocation, and vaccine data when available.
```mermaid
flowchart LR
subgraph "Final Full Compiled Data Sets"
  H[Lauring Lab Influenza A/B Sequencing, Geolocation, Vaccine, and Demographic Data]:::final1
  I[Lauring Lab Influenza A/B Sequencing that only has matches to Demographic data and then any vaccine and geolocation data that matches this filtered data set]:::final2
  classDef final1 stroke:#0f0,stroke-width:4px
  classDef final2 stroke:#00f,stroke-width:4px
end
subgraph "Subset Data Sets"
  E[Lauring Lab Influenza A/B Sequencing Data and Geolocation Data]
  F[Lauring Lab Influenza A/B Sequencing Data and Vaccine Data]
  G[Lauring Lab Influenza A/B Sequencing Data and Demographic Data]
end
subgraph "Data Inputs"
  A[Lauring Lab Influenza A/B Sequencing Data]
  B[Vaccine Data from Lighthouse]
  C[Geolocation Data from ArcGIS Pro 3.6]
  D[Demographic Data from partnering Stakeholders]
end
 %% to add a comment to this code I need to add these 2 % and then anything after that is a comment until a new line %%

```
