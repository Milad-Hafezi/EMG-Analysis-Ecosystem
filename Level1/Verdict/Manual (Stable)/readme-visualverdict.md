# Visual-Verdict

## Description
Visual-Verdict is a spreadsheet-based tool for systematically recording human assessments of EMG responses to ICMS. It provides a structured framework for researchers to evaluate EMG activity observed in the PDF reports generated by the Raw Signal Plotter, determining threshold currents for each stimulation site and muscle combination.

## Inputs
- PDF reports generated by the Raw Signal Plotter
- Human visual assessment

## Outputs
- Excel spreadsheet (.xlsx) with 17 sheets:
  - 16 sheets (one per muscle) containing threshold assessments
    - Column D records threshold current values for each stimulation site
    - Column F provides space for comments and observations
  - 1 metadata sheet (Sheet 17) containing:
    - Case-specific information (case number in cell F1)
    - Muscle abbreviations and definitions (using lookup table in columns J and K)
    - Channel number assignments (column A)
    - Additional experiment parameters

## Usage
1. Make a copy of the Sample_Verdict template spreadsheet
2. Review the PDF reports from the Raw Signal Plotter
3. For each muscle sheet (1-16), fill in threshold values in Column D for each stimulation site
4. Add any comments or observations in Column F
5. Complete Sheet 17 with:
   - Case number in cell F1
   - Muscle abbreviations using the lookup table
   - Channel number assignments in column A
6. Save the spreadsheet with appropriate naming for version tracking

## Dependencies
- Microsoft Excel or compatible spreadsheet software
- PDF viewer for reviewing Raw Signal Plotter reports

## Notes
- Always work from a copy of the template to preserve the original format
- Maintain consistent naming conventions for version tracking
- The completed spreadsheet serves as critical input for the MetaData Generator
- Threshold values documented here will determine which responses are included in topographic mapping
- Comments can provide valuable context for unusual responses or potential artifacts
