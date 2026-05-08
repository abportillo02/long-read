from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import inch

output_path = "/mnt/user-data/outputs/KZFP_AMPure_Size_Selection_Protocol.pdf"
doc = SimpleDocTemplate(output_path, pagesize=letter,
                        leftMargin=0.75*inch, rightMargin=0.75*inch,
                        topMargin=0.75*inch, bottomMargin=0.75*inch)

styles = getSampleStyleSheet()

# Custom styles
title_style = ParagraphStyle('CustomTitle', parent=styles['Title'],
    fontSize=16, spaceAfter=4, textColor=colors.HexColor('#003366'))

subtitle_style = ParagraphStyle('Subtitle', parent=styles['Normal'],
    fontSize=10, spaceAfter=12, textColor=colors.HexColor('#555555'), italic=True)

heading_style = ParagraphStyle('SectionHeading', parent=styles['Heading2'],
    fontSize=12, spaceBefore=14, spaceAfter=4,
    textColor=colors.white,
    backColor=colors.HexColor('#003366'),
    leftPadding=6, borderPad=4)

subheading_style = ParagraphStyle('SubHeading', parent=styles['Heading3'],
    fontSize=10, spaceBefore=8, spaceAfter=3,
    textColor=colors.HexColor('#003366'))

body_style = ParagraphStyle('Body', parent=styles['Normal'],
    fontSize=9.5, spaceAfter=4, leading=14)

step_style = ParagraphStyle('Step', parent=styles['Normal'],
    fontSize=9.5, spaceAfter=5, leading=14, leftIndent=18)

warning_style = ParagraphStyle('Warning', parent=styles['Normal'],
    fontSize=9, spaceAfter=6, leading=13, leftIndent=12,
    textColor=colors.HexColor('#7B0000'),
    backColor=colors.HexColor('#FFF3F3'),
    borderColor=colors.HexColor('#CC0000'),
    borderWidth=1, borderPad=5)

tip_style = ParagraphStyle('Tip', parent=styles['Normal'],
    fontSize=9, spaceAfter=6, leading=13, leftIndent=12,
    textColor=colors.HexColor('#004400'),
    backColor=colors.HexColor('#F0FFF0'),
    borderColor=colors.HexColor('#007700'),
    borderWidth=1, borderPad=5)

note_style = ParagraphStyle('Note', parent=styles['Normal'],
    fontSize=9, spaceAfter=6, leading=13, leftIndent=12,
    textColor=colors.HexColor('#1a1a5e'),
    backColor=colors.HexColor('#F0F4FF'),
    borderColor=colors.HexColor('#3355aa'),
    borderWidth=1, borderPad=5)

story = []

# Title
story.append(Paragraph("AMPure XP Size Selection Protocol", title_style))
story.append(Paragraph("Single-sided SPRI cleanup — Remove fragments &lt;~1,000 bp | Adjusted for 73 µL input", subtitle_style))
story.append(Paragraph("Application: Post-capture PCR cleanup for KZFP gene hybridization capture", subtitle_style))
story.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor('#003366'), spaceAfter=10))

# Purpose box
story.append(Paragraph(
    "<b>Purpose:</b> Remove short fragments (&lt;~1,000 bp) from post-capture long-range PCR product "
    "prior to third-generation sequencing. Uses a 0.45x AMPure XP bead ratio to retain fragments in "
    "the KZFP gene size range (~1,200–15,800 bp) while discarding primer dimers, adapters, and small fragments.",
    note_style))
story.append(Spacer(1, 6))

# Materials
story.append(Paragraph("Materials Required", heading_style))
story.append(Spacer(1, 4))

mat_data = [
    ["Reagent / Equipment", "Specification", "Notes"],
    ["Agencourt AMPure XP Beads", "Standard (Beckman Coulter A63880)", "Equilibrate to RT before use"],
    ["Ethanol", "200-proof/absolute, freshly diluted to 80% v/v", "Prepare fresh on day of use"],
    ["Elution Buffer", "10 mM Tris-HCl, pH 8.5 (e.g., Qiagen Buffer EB)", "Equilibrate to RT before use"],
    ["Nuclease-free water", "—", "For bead ratio calculation if needed"],
    ["Low-bind 1.5 mL tubes", "—", "Use throughout"],
    ["Magnetic rack", "1.5 mL tube compatible", "—"],
    ["Microcentrifuge", "—", "For pulse-spins"],
    ["Pipettes", "Ranging 1–1000 µL, with aerosol-resistant tips", "Use low-retention tips if available"],
    ["TapeStation / Bioanalyzer", "gDNA ScreenTape or equivalent", "For QC after cleanup"],
]

mat_table = Table(mat_data, colWidths=[2.1*inch, 2.4*inch, 2.1*inch])
mat_table.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#003366')),
    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 8.5),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.HexColor('#F5F8FF'), colors.white]),
    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#AAAACC')),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 4),
    ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ('LEFTPADDING', (0,0), (-1,-1), 5),
]))
story.append(mat_table)
story.append(Spacer(1, 8))

# Volume summary table
story.append(Paragraph("Volume Summary", heading_style))
story.append(Spacer(1, 4))

vol_data = [
    ["Parameter", "Value", "Basis"],
    ["Input sample volume", "73 µL", "Your post-capture PCR product"],
    ["AMPure XP beads (0.45x)", "33 µL", "73 µL × 0.45 = 32.85 µL, rounded to 33 µL"],
    ["Size cutoff (lower)", "~1,000–1,200 bp", "Fragments below this are discarded"],
    ["Ethanol wash", "1 mL × 2 washes", "80% v/v, freshly prepared"],
    ["Elution buffer", "32 µL", "Scaled down to maximize concentration"],
    ["Expected recovery volume", "30 µL", "Transfer supernatant after elution"],
]

vol_table = Table(vol_data, colWidths=[2.2*inch, 1.5*inch, 2.9*inch])
vol_table.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#003366')),
    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTNAME', (0,1), (0,-1), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 8.5),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.HexColor('#F5F8FF'), colors.white]),
    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#AAAACC')),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 4),
    ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ('LEFTPADDING', (0,0), (-1,-1), 5),
]))
story.append(vol_table)
story.append(Spacer(1, 8))

# Important notes box
story.append(Paragraph(
    "<b>⚠ Important — Read before starting:</b><br/>"
    "• Equilibrate AMPure XP beads AND Elution Buffer to room temperature for at least 30 minutes before use<br/>"
    "• Vortex beads until completely homogeneous immediately before pipetting<br/>"
    "• All steps performed at room temperature unless otherwise stated<br/>"
    "• Use low-bind tubes throughout to minimize DNA loss<br/>"
    "• Do NOT use the custom size-selection beads from the IDT reagent prep — use standard AMPure XP here<br/>"
    "• Prepare fresh 80% ethanol on the day of use",
    warning_style))
story.append(Spacer(1, 6))

# Protocol steps
story.append(Paragraph("Protocol", heading_style))
story.append(Spacer(1, 4))

# Step 1
story.append(Paragraph("Step 1 — Bead Binding", subheading_style))
story.append(Paragraph(
    "<b>1.</b> Transfer <b>73 µL</b> of post-capture PCR product into a new <b>1.5 mL low-bind tube</b>.",
    step_style))
story.append(Paragraph(
    "<b>2.</b> Add <b>33 µL</b> of AMPure XP beads (0.45x ratio) directly to the sample.",
    step_style))
story.append(Paragraph(
    "<b>3.</b> Mix thoroughly by carefully pipetting up and down <b>10 times</b> until the solution is completely homogeneous. "
    "Avoid introducing air bubbles. Use a 100 µL pipette tip — do not use tips smaller than 100 µL as shear forces can damage long fragments.",
    step_style))
story.append(Paragraph(
    "<b>4.</b> Incubate for <b>7 minutes</b> at room temperature. Gently flick or tap the tube every 1–2 minutes to keep beads in suspension.",
    step_style))
story.append(Paragraph(
    "<b>5.</b> Pulse-spin the tube briefly in a microcentrifuge to collect liquid.",
    step_style))
story.append(Paragraph(
    "<b>6.</b> Place tube on the magnetic rack and wait until the solution appears completely clear (~2 minutes).",
    step_style))
story.append(Paragraph(
    "<b>7.</b> Carefully remove and <b>discard the supernatant</b> without disturbing the bead pellet. "
    "The supernatant contains fragments &lt;~1,000 bp — this is waste.",
    step_style))

story.append(Spacer(1, 4))
story.append(Paragraph("Step 2 — Ethanol Washes (×2)", subheading_style))
story.append(Paragraph(
    "<b>8.</b> Keep the tube on the magnetic rack. Slowly dispense <b>1 mL of freshly prepared 80% ethanol</b> "
    "against the tube wall <i>opposite</i> the bead pellet. <b>Do not resuspend the beads.</b>",
    step_style))
story.append(Paragraph(
    "<b>9.</b> Wait <b>30 seconds</b>, then carefully remove and discard the ethanol.",
    step_style))
story.append(Paragraph(
    "<b>10.</b> Repeat steps 8–9 once more for a total of <b>2 ethanol washes</b>.",
    step_style))
story.append(Paragraph(
    "<b>11.</b> Pulse-spin the tube, place back on the magnet, and remove any residual ethanol from the bottom using a fresh pipette tip.",
    step_style))
story.append(Paragraph(
    "<b>12.</b> Check for any remaining droplets on the tube walls. If present, pulse-spin again and remove.",
    step_style))

story.append(Spacer(1, 4))
story.append(Paragraph("Step 3 — Bead Drying", subheading_style))
story.append(Paragraph(
    "<b>13.</b> Leave the tube on the magnetic rack with the cap open and air-dry the beads for <b>60–120 seconds</b>.",
    step_style))
story.append(Paragraph(
    "<b>Checkpoint:</b> The bead pellet is ready when it just begins to look matte rather than shiny/wet — "
    "add elution buffer immediately at this point. Do NOT wait until the pellet cracks.",
    tip_style))

story.append(Spacer(1, 4))
story.append(Paragraph("Step 4 — Elution", subheading_style))
story.append(Paragraph(
    "<b>14.</b> Remove the tube from the magnetic rack. Add <b>32 µL of Elution Buffer</b> (10 mM Tris-HCl, pH 8.5) "
    "directly onto the bead pellet.",
    step_style))
story.append(Paragraph(
    "<b>15.</b> Pipette slowly up and down until beads are <b>completely and uniformly resuspended</b>. "
    "Take your time — thorough resuspension is critical for maximum yield.",
    step_style))
story.append(Paragraph(
    "<b>16.</b> Incubate at <b>37°C for 5 minutes</b> to maximize elution yield. "
    "(Room temperature incubation for 5 minutes is also acceptable but 37°C is recommended given the smaller input volume.)",
    step_style))
story.append(Paragraph(
    "<b>17.</b> Pulse-spin the tube briefly, then place back on the magnetic rack.",
    step_style))
story.append(Paragraph(
    "<b>18.</b> Wait approximately <b>2 minutes</b> until the solution appears completely clear.",
    step_style))
story.append(Paragraph(
    "<b>19.</b> Carefully transfer <b>30 µL of the supernatant</b> to a new 1.5 mL low-bind tube. "
    "Avoid disturbing the bead pellet.",
    step_style))
story.append(Paragraph(
    "<b>20.</b> Discard the beads.",
    step_style))

story.append(Spacer(1, 4))
story.append(Paragraph("Step 5 — Quality Control", subheading_style))
story.append(Paragraph(
    "<b>21.</b> Quantify your DNA using a <b>Qubit fluorometer</b> (dsDNA HS assay) or equivalent fluorometric method.",
    step_style))
story.append(Paragraph(
    "<b>22.</b> Verify fragment size distribution using <b>TapeStation with gDNA ScreenTape</b> or equivalent. "
    "Confirm that fragments below ~1,000 bp have been removed and your KZFP-length target fragments are intact.",
    step_style))

story.append(Spacer(1, 6))
story.append(Paragraph(
    "<b>💡 Tip — Over-drying warning:</b> Since your input volume is smaller than the standard IDT protocol, "
    "your bead pellet will be relatively small. Be especially vigilant not to over-dry in step 13 — "
    "over-drying long fragments makes them very difficult to elute and can significantly reduce yield. "
    "When in doubt, add elution buffer sooner rather than later.",
    tip_style))

story.append(Spacer(1, 6))

# Storage
story.append(Paragraph("Storage", heading_style))
story.append(Spacer(1, 4))
stor_data = [
    ["Timeframe", "Storage condition"],
    ["Up to 24 hours", "4°C"],
    ["Up to 4 weeks", "-20°C"],
    ["Long term", "-80°C (recommended for precious samples)"],
]
stor_table = Table(stor_data, colWidths=[2.5*inch, 4.1*inch])
stor_table.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#003366')),
    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 9),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.HexColor('#F5F8FF'), colors.white]),
    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#AAAACC')),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 5),
    ('BOTTOMPADDING', (0,0), (-1,-1), 5),
    ('LEFTPADDING', (0,0), (-1,-1), 6),
]))
story.append(stor_table)
story.append(Spacer(1, 8))

# Troubleshooting
story.append(Paragraph("Troubleshooting", heading_style))
story.append(Spacer(1, 4))
ts_data = [
    ["Issue", "Possible Cause", "Solution"],
    ["Low yield after elution", "Beads over-dried", "Reduce air-dry time to 60 sec; elute at 37°C"],
    ["Small fragments still present", "Bead ratio too high (>0.45x)", "Repeat cleanup with fresh 0.45x beads"],
    ["Large fragments lost", "Bead ratio too low (<0.45x)", "Check pipetting accuracy; re-calibrate"],
    ["Poor bead resuspension", "Ethanol carryover / over-drying", "Ensure all ethanol removed; resuspend thoroughly"],
    ["Very low DNA concentration", "Small input volume + losses", "Reduce elution volume to 20 µL to concentrate"],
]
ts_table = Table(ts_data, colWidths=[1.8*inch, 2.1*inch, 2.7*inch])
ts_table.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#003366')),
    ('TEXTCOLOR', (0,0), (-1,0), colors.white),
    ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
    ('FONTSIZE', (0,0), (-1,-1), 8.5),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.HexColor('#F5F8FF'), colors.white]),
    ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#AAAACC')),
    ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ('TOPPADDING', (0,0), (-1,-1), 4),
    ('BOTTOMPADDING', (0,0), (-1,-1), 4),
    ('LEFTPADDING', (0,0), (-1,-1), 5),
    ('WORDWRAP', (0,0), (-1,-1), True),
]))
story.append(ts_table)
story.append(Spacer(1, 8))

# Footer note
story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor('#AAAACC'), spaceAfter=6))
story.append(Paragraph(
    "Adapted from: IDT High-throughput Hybridization Capture Enrichment of Long Genomic Fragments Protocol (RA-DOC-614 REV 02), "
    "Section 9B Clean-Up. Bead ratio and volumes adjusted for 73 µL post-capture PCR input. "
    "For Research Use Only.",
    ParagraphStyle('Footer', parent=styles['Normal'], fontSize=7.5, textColor=colors.grey)))

doc.build(story)
print("PDF created successfully.")