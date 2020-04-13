// The TStyle
extern TStyle t2kStyle;

// The rootlogon
void rootlogon() {

  t2kStyle = new TStyle("t2kStyle","T2K approved plots style");

  // Plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetFrameLineWidth(2);

  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetCanvasBorderSize(0);    
  t2kStyle->SetCanvasColor(0);

  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadBorderSize(0);    
  t2kStyle->SetPadColor(0);

  t2kStyle->SetTitleFillColor(0);

  t2kStyle->SetFillStyle(0);

  t2kStyle->SetStatColor(0);

  // Grids
  t2kStyle->SetPadGridX(1);
  t2kStyle->SetPadGridY(1);
  t2kStyle->SetGridStyle(3);
  t2kStyle->SetGridColor(0);
  t2kStyle->SetGridWidth(1);

  // Legends
  t2kStyle->SetLegendBorderSize(0); 
  t2kStyle->SetLegendFillColor(0);
  t2kStyle->SetLegendFont(42);
  //t2kStyle->SetLegendTextSize(5);


  // Paper and margin
  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.10);
  t2kStyle->SetPadRightMargin(0.1);
  t2kStyle->SetPadBottomMargin(0.16);
  t2kStyle->SetPadLeftMargin(0.1);

  // Large Helvetica fonts
  t2kStyle->SetTextFont(42);
  t2kStyle->SetTextSize(0.08);

  // Axes
  t2kStyle->SetLabelFont(42,"xyz");
  t2kStyle->SetLabelSize(0.04,"xyz");
  t2kStyle->SetLabelOffset(0.015, "xyz");
  t2kStyle->SetLabelColor(1, "xyz");

  t2kStyle->SetAxisColor(1);

  // Titles
  t2kStyle->SetTitleFont(42,"xyz");
  t2kStyle->SetTitleSize(0.05,"xyz");
  t2kStyle->SetTitleOffset(1.1, "yz");
  t2kStyle->SetTitleOffset(1.0, "x");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleColor(1, "xyzt");
  t2kStyle->SetTitleX(0.50);
  t2kStyle->SetTitleY(1.00);

  // Align title with centre
  t2kStyle->SetTitleAlign(23);
  t2kStyle->SetTitleBorderSize(0);    

  // Stats box
  t2kStyle->SetStatFont(42);
  //t2kStyle->SetStatFontSize(0.07);
  t2kStyle->SetStatBorderSize(0);


  // Default Histogram properties
  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetMarkerSize(1.0);
  t2kStyle->SetHistLineWidth(2);

  //t2kStyle->SetHistLineColor(kBlack);
  //t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // Default Function properties
  //t2kStyle->SetFuncStyle(1);
  //t2kStyle->SetFuncColor(kRed);
  t2kStyle->SetFuncWidth(2);

  // Change Error formatting
  //t2kStyle->SetErrorX(0.5);

  // SetOpt for all
  //t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);
  t2kStyle->SetOptDate(0);
  t2kStyle->SetDateX(0.01);
  t2kStyle->SetDateY(0.01);

  // Tick marks on top and RHS of plots
  t2kStyle->SetTickLength(0.03, "xyzt");
  t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPaintTextFormat("1.2f");

  SetPalette();

  gROOT->SetStyle("t2kStyle");

  // gROOT->ForceStyle();
}

// Make the TPalette
static void SetPalette() {

  /*
  // Add a greyscale palette for 2D plots
  int ncol=50;
  double dcol = 1./float(ncol);
  double gray = 1;

  TColor **theCols = new TColor*[ncol];

  for (int i=0;i<ncol;i++) theCols[i] = new TColor(999-i,0.0,0.7,0.7);

  for (int j = 0; j < ncol; j++) {
  theCols[j]->SetRGB(gray,gray,gray);
  gray -= dcol;
  }

  int ColJul[100];
  for  (int i=0; i<100; i++) ColJul[i]=999-i;
  t2kStyle->SetPalette(ncol,ColJul);
   */

  // Define a nicer color palette (red->blue)

  // Uncomment these lines for a color palette (default is B&W)
  //    t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  //    const Int_t NRGBs = 5;
  //    const Int_t NCont = 255;

  //    // Uncomment these colours for the rainbow (blue -> yellow -> red) palette
  //    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };

  //    // Uncomment these to use the black-bady palette
  //    Double_t white[3]   = {1, 1, 1};
  //    Double_t yellow[3]  = {1, 1, 0};
  //    Double_t orange[3]  = {1, 0.5, 0};
  //    Double_t crimson[3]     = {1, 0, 0};
  //    Double_t black[3]   = {0, 0, 0};
  //    Double_t red[NRGBs]    = {black[0], crimson[0], orange[0], yellow[0], white[0]};
  //    Double_t green[NRGBs]  = {black[1], crimson[1], orange[1], yellow[1], white[1]};
  //    Double_t blue[NRGBs]   = {black[2], crimson[2], orange[2], yellow[2], white[2]};
  //    Double_t stops[NRGBs]  = {   0.00,    0.25,       0.5,      0.75,     1.00};

  // Uncomment these to use the blue->white->red palette (good for correlation matrices)
  TColor::InitializeColors();
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  /*
     way this works is that you define the stops and what fraction of colour
     should be applied. e.g. for stop 0.10 we have 0.25 red, 0.25 green, 1.00
     blue
   */

  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };

  int start = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  t2kStyle->SetNumberContours(NCont);

}
