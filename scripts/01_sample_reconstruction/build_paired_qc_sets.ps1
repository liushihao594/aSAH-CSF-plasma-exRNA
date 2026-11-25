param(
  [Parameter(Mandatory=$true)][string]$ValidatedLibraryQcCsv,
  [Parameter(Mandatory=$true)][string]$PrimaryPairsCsv,
  [Parameter(Mandatory=$true)][string]$OutputDirectory
)

# Public, path-parameterized release copy of the Stage 0C pairing/QC logic.
# It writes patient-level local outputs. Keep OutputDirectory under an ignored
# data_private/ or derived_local/ path and do not commit those outputs.
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
$qc = Import-Csv -LiteralPath $ValidatedLibraryQcCsv
$pairs = Import-Csv -LiteralPath $PrimaryPairsCsv
$mapped = @{}
foreach ($r in $qc) { $mapped[$r.subject + '::' + $r.label + '::' + $r.biofluid] = $r }

$result = foreach ($p in $pairs) {
  $csf = $mapped[$p.subject + '::' + $p.collection_label + '::CSF']
  $plasma = $mapped[$p.subject + '::' + $p.collection_label + '::plasma']
  if (-not $csf -or -not $plasma) { throw 'Missing member of original pair' }
  [pscustomobject]@{
    pair_id = $p.subject + '::' + $p.collection_label
    subject = $p.subject
    collection_label = $p.collection_label
    centered_label = [int]$p.centered_label
    CSF_GSM = $csf.GSM
    plasma_GSM = $plasma.GSM
    CSF_expression_id = $csf.expression_id
    plasma_expression_id = $plasma.expression_id
    CSF_QC = $csf.QC
    plasma_QC = $plasma.QC
    CSF_transcriptome_reads = $csf.transcriptome_reads
    plasma_transcriptome_reads = $plasma.transcriptome_reads
    CSF_reference_genome_reads = $csf.reference_genome_reads
    plasma_reference_genome_reads = $plasma.reference_genome_reads
    both_QC_PASS = ($csf.QC -eq 'PASS' -and $plasma.QC -eq 'PASS')
  }
}

if ($result.Count -ne 153) { throw "Expected 153 complete pairs; observed $($result.Count)" }
if (@($result | Where-Object {$_.both_QC_PASS}).Count -ne 141) { throw 'Expected 141 QC PASS-only pairs' }

$result | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'FULL_153_pairs.csv') -NoTypeInformation -Encoding utf8
$result | Where-Object {$_.both_QC_PASS} | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'QC_PASS_ONLY_141_pairs.csv') -NoTypeInformation -Encoding utf8
$result | Where-Object {-not $_.both_QC_PASS} | Export-Csv -LiteralPath (Join-Path $OutputDirectory 'QC_AFFECTED_PAIRS.csv') -NoTypeInformation -Encoding utf8
Write-Output "PASS: $($result.Count) complete pairs; $(@($result | Where-Object {$_.both_QC_PASS}).Count) QC PASS-only pairs"
