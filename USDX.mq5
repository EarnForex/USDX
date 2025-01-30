//+------------------------------------------------------------------+
//|                                                   USDX Indicator |
//|                                 Copyright © 2012-2025, EarnForex |
//|                                        https://www.earnforex.com |
//+------------------------------------------------------------------+
#property copyright "www.EarnForex.com, 2012-2025"
#property link      "https://www.earnforex.com/metatrader-indicators/USDX/"
#property version   "1.04"

#property description "USDX - indicator of the US Dollar Index."
#property description "Displays a USDX chart in a separate window of the current chart."
#property description "Based on EUR/USD, USD/JPY, GBP/USD, USD/CAD, USD/CHF and USD/SEK."
#property description "All these pairs should be added to Market Watch for indicator to work."
#property description "Two customizable moving averages can be applied to the index."
#property description "Can be easily modified via input parameters to calculate any currency index."

#property indicator_separate_window

#property indicator_plots 3
#property indicator_buffers 7

#property indicator_type1  DRAW_COLOR_CANDLES
#property indicator_color1 clrDodgerBlue, clrRed
#property indicator_label1 "USDX"

#property indicator_type2  DRAW_LINE
#property indicator_color2 clrWhite
#property indicator_label2 "MA1"

#property indicator_type3  DRAW_LINE
#property indicator_color3 clrLimeGreen
#property indicator_label3 "MA2"

input int MaxBars = 500; // MaxBars (Set to 0 to disable)
input string IndexPairs = "EURUSD, USDJPY, GBPUSD, USDCAD, USDSEK, USDCHF"; // IndexPairs (Update only currency pairs' names if they are different from your broker's)
input string IndexCoefficients = "-0.576, 0.136, -0.119, 0.091, 0.042, 0.036"; // IndexCoefficients (Do not change for USDX)
input double IndexInitialValue = 50.14348112; // IndexInitialValue (Do not change for USDX)
input int MA_Period1 = 13; // MA_Period1 (Set to 0 to disable this MA)
input int MA_Period2 = 17; // MA_Period2 (Set to 0 to disable this MA)
input ENUM_MA_METHOD MA_Mode1 = MODE_SMA;
input ENUM_MA_METHOD MA_Mode2 = MODE_SMA;
input ENUM_APPLIED_PRICE MA_PriceType1 = PRICE_CLOSE;
input ENUM_APPLIED_PRICE MA_PriceType2 = PRICE_CLOSE;

// Data structure for pairs' data
class CPairData
{
public:
    double   Open[], High[], Low[], Close[];
    datetime Time[];
};

// Indicator buffers:
double Open[], High[], Low[], Close[];
double Color[];
double MA1[], MA2[];

// Pairs' arrays:
string Pairs[];
double Coefficients[];
bool FirstRun[];
CPairData PairData[];
int PairsLastBars[];
int Pairs_k[];
int Pairs_NotFound[]; // 0 - OK, 1 - first failure, 2 - final failure.

int LastBars = 0, RatesTotal = 0, LastMA = 0;
bool AllFound = true;

// Gets pair names and coefficients out of the input parameters and creates the respective arrays.
int InitializePairs()
{
    int n = 0, count = 0, coef_n = 0, coef_count = 0;
    // Count elements.
    while (n != -1)
    {
        n = StringFind(IndexPairs, ",", n + 1);
        coef_n = StringFind(IndexCoefficients, ",", coef_n + 1);
        if (n > -1) count++;
        if (coef_n > -1) coef_count++;
    }
    count++;
    coef_count++;

    if (count != coef_count)
    {
        Alert("Each currency pair of the Index should have a corresponding coefficient.");
        return -1;
    }

    ArrayResize(Pairs, count);
    ArrayResize(PairsLastBars, count);
    ArrayResize(Pairs_k, count);
    ArrayResize(Pairs_NotFound, count);
    ArrayResize(Coefficients, count);
    ArrayResize(PairData, count);
    ArrayResize(FirstRun, count);
    ArrayInitialize(FirstRun, true);
    ArrayInitialize(PairsLastBars, 0);
    ArrayInitialize(Pairs_NotFound, 0);

    n = 0;
    coef_n = 0;
    int prev_n = 0, prev_coef_n = 0;
    for (int i = 0; i < count; i++)
    {
        n = StringFind(IndexPairs, ",", n + 1);
        coef_n = StringFind(IndexCoefficients, ",", coef_n + 1);
        // To remove the trailing comma
        if (prev_n > 0) prev_n++;
        if (prev_coef_n > 0) prev_coef_n++;
        Pairs[i] = StringSubstr(IndexPairs, prev_n, n - prev_n);
        Coefficients[i] = StringToDouble(StringSubstr(IndexCoefficients, prev_coef_n, coef_n - prev_coef_n));
        StringTrimLeft(Pairs[i]);
        StringTrimRight(Pairs[i]);
        Print(Pairs[i], ": ", Coefficients[i]);
        if (!(bool)SymbolInfoInteger(Pairs[i], SYMBOL_SELECT)) SymbolSelect(Pairs[i], true);
        if (!SymbolInfoDouble(Pairs[i], SYMBOL_BID))
        {
            Print(Pairs[i], " not found. Could be a delay while selecting it to Market Watch. Will try again inside OnCalculate().");
            AllFound = false;
        }
        else
        {
            Print(Pairs[i], " found.");
        }
        prev_n = n;
        prev_coef_n = coef_n;
    }

    return 0;
}

int OnInit()
{
    if (StringLen(IndexPairs) == 0)
    {
        Alert("Please enter the USDX currency pairs (EUR/USD, USD/JPY, GBP/USD, USD/CAD, USD/CHF and USD/SEK) into IndexPairs input parameter. Use the currency pairs' names as they appear in your market watch.");
    }

    if (StringLen(IndexCoefficients) == 0)
    {
        Alert("Please enter the the Index coefficients as coma-separated values in IndexCoefficients input parameter.");
    }

    if (InitializePairs() < 0)
    {
        Print("Couldn't initialize pairs.");
    }

    SetIndexBuffer(0, Open, INDICATOR_DATA);
    SetIndexBuffer(1, High, INDICATOR_DATA);
    SetIndexBuffer(2, Low, INDICATOR_DATA);
    SetIndexBuffer(3, Close, INDICATOR_DATA);
    SetIndexBuffer(4, Color, INDICATOR_COLOR_INDEX);
    if (MA_Period1 > 0)
    {
        SetIndexBuffer(5, MA1, INDICATOR_DATA);
        PlotIndexSetDouble(1, PLOT_EMPTY_VALUE, EMPTY_VALUE);
        PlotIndexSetString(1, PLOT_LABEL, "MA(" + IntegerToString(MA_Period1) + ")");
        PlotIndexSetInteger(1, PLOT_DRAW_TYPE, DRAW_LINE);
    }
    else PlotIndexSetInteger(1, PLOT_DRAW_TYPE, DRAW_NONE);
    if (MA_Period2 > 0)
    {
        SetIndexBuffer(6, MA2, INDICATOR_DATA);
        PlotIndexSetDouble(2, PLOT_EMPTY_VALUE, EMPTY_VALUE);
        PlotIndexSetString(2, PLOT_LABEL, "MA(" + IntegerToString(MA_Period2) + ")");
        PlotIndexSetInteger(2, PLOT_DRAW_TYPE, DRAW_LINE);
    } else PlotIndexSetInteger(2, PLOT_DRAW_TYPE, DRAW_NONE);
    
    IndicatorSetInteger(INDICATOR_DIGITS, 4);
    IndicatorSetString(INDICATOR_SHORTNAME, "USDX");

    PlotIndexSetDouble(0, PLOT_EMPTY_VALUE, EMPTY_VALUE);

    return INIT_SUCCEEDED;
}

int OnCalculate(const int rates_total,
                const int prev_calculated,
                const datetime &Time[],
                const double &open[],
                const double &high[],
                const double &low[],
                const double &close[],
                const long &tick_volume[],
                const long &volume[],
                const int &spread[])
{
    if (!AllFound) // Making sure that the Terminal has all the required currency pairs.
    {
        AllFound = true;
        for (int i = 0; i < ArraySize(Pairs); i++)
        {
            if (!SymbolInfoDouble(Pairs[i], SYMBOL_BID))
            {
                AllFound = false;
                if (Pairs_NotFound[i] == 2) return 0; // Cannot calculate if one of the pairs is missing.
                Pairs_NotFound[i]++;
                if (Pairs_NotFound[i] == 1) Print(Pairs[i], " not found. Still waiting due to a potential delay while selecting it to Market Watch.");
                else if (Pairs_NotFound[i] == 2) Alert(Pairs[i], " not found. Please check its naming and whether your broker offers such a pair.");
            }
            if ((Pairs_NotFound[i] == 1) && (AllFound)) // If this pair was missing at some point, but is now OK.
            {
                Print(Pairs[i], " found.");
            }
        }
    }
    RatesTotal = rates_total;

    // Clear all data for the first run or when one of the currency pairs fails to load.
    if (prev_calculated == 0)
    {
        for (int i = 0; i < rates_total; i++)
        {
            Open[i] = EMPTY_VALUE;
            High[i] = EMPTY_VALUE;
            Low[i] = EMPTY_VALUE;
            Close[i] = EMPTY_VALUE;
            if (MA_Period1 > 0) MA1[i] = EMPTY_VALUE;
            if (MA_Period2 > 0) MA2[i] = EMPTY_VALUE;
        }
    }

    int pairs_count = ArraySize(Pairs);
    int NewLastBars = -1;

    for (int i = 0; i < pairs_count; i++)
    {
        UpdateArrays(i);
        if (PairsLastBars[i] > 0) FirstRun[i] = false; // Initial array filling complete.

        // Charts not ready yet.
        if ((ArraySize(PairData[i].Time) == 0) || (ArraySize(PairData[i].Open) == 0) || (ArraySize(PairData[i].High) == 0) || (ArraySize(PairData[i].Low) == 0) || (ArraySize(PairData[i].Close) == 0))
        {
            Print(Pairs[i], " failed to load. T:", ArraySize(PairData[i].Time), " O:", ArraySize(PairData[i].Open), " H:", ArraySize(PairData[i].High), " L:", ArraySize(PairData[i].Low), " C:", ArraySize(PairData[i].Close));
            // Reset for all pairs because everything has to be re-calculated.
            ArrayInitialize(PairsLastBars, 0);
            LastBars = 0;
            return prev_calculated;
        }

        // Data arrays are filled.

        // Starting point for iterations through arrays.
        Pairs_k[i] = ArraySize(PairData[i].Time) - 1;

        // Skip currency pair if its arrays are empty for some reason.
        if (Pairs_k[i] < 0) continue;

        // Go through OHLC arrays syncing them by Time and calculate indexes.
        // Leave EMPTY_VALUE for desynced bars.
        for (int j = RatesTotal - 1; j >= LastBars; j--)
        {
            // Got to the end of currency pairs arrays but still have new bars in the current chart.
            if (Pairs_k[i] < 0)
            {
                // Assign EV to the new bars without equivalent on index currency pairs
                Open[j] = EMPTY_VALUE;
                High[j] = EMPTY_VALUE;
                Low[j] = EMPTY_VALUE;
                Close[j] = EMPTY_VALUE;
                continue;
            }
            // Initialize indicator with IndexInitialValue (50.14348112 for USDX) on first currency pair (to multiply everything with it)
            if (i == 0)
            {
                Open[j] = IndexInitialValue;
                High[j] = IndexInitialValue;
                Low[j] = IndexInitialValue;
                Close[j] = IndexInitialValue;
            }
            // On second and further currency pairs, have to check if indicator value was already set to EV, if so - no need to check this i.
            else if (Open[j] == EMPTY_VALUE) continue;
            // Cycle through OHLC arrays of the particular currency pair
            for (int n = Pairs_k[i]; n >= 0; n--)
            {
                // Got to the right bar.
                if (PairData[i].Time[n] == Time[j])
                {
                    Open[j] *= MathPow(PairData[i].Open[n], Coefficients[i]);
                    High[j] *= MathPow(PairData[i].High[n], Coefficients[i]);
                    Low[j] *= MathPow(PairData[i].Low[n], Coefficients[i]);
                    Close[j] *= MathPow(PairData[i].Close[n], Coefficients[i]);
                    // Decrease counter only if the Time value was found, otherwise the current element can be the next one sought.
                    Pairs_k[i] = n - 1;
                    // Last pair.
                    if (i == pairs_count - 1)
                    {
                        // Once all pairs are in, assign a color.
                        if (Close[j] >= Open[j]) Color[j] = 0;
                        else Color[j] = 1;
                        // LastBar needs to be updated if it's latest bar of the 6th pair.
                        // Do it via NewLastBars for not to mess with the cycle.
                        if (j == RatesTotal - 1) NewLastBars = RatesTotal - 1;

                        // Reached last bar? Calculate MAs!
                        if (j == LastBars) CalculateMAs();
                    }
                    break;
                }
                // Got too far back in time, the bar is missing for this currency pair.
                else if (PairData[i].Time[n] < Time[j])
                {
                    // Try using outdated data while the new bar is arriving for the pair:
                    Open[j] *= MathPow(PairData[i].Open[n], Coefficients[i]);
                    High[j] *= MathPow(PairData[i].High[n], Coefficients[i]);
                    Low[j] *= MathPow(PairData[i].Low[n], Coefficients[i]);
                    Close[j] *= MathPow(PairData[i].Close[n], Coefficients[i]);
                    break;
                }
                // If got to the last element without finding, the bar is missing.
                else if (n == 0)
                {
                    Open[j] = EMPTY_VALUE;
                    High[j] = EMPTY_VALUE;
                    Low[j] = EMPTY_VALUE;
                    Close[j] = EMPTY_VALUE;
                }
            }
        }
        // Finally updating LastBars if got to the last bar of the last currency pair inside the cycle.
        if (NewLastBars > -1) LastBars = NewLastBars;
    }

    // If limit for max number of bars is set, apply it.
    if ((MaxBars > 0) && (LastBars - RatesTotal > MaxBars)) LastBars = RatesTotal - MaxBars;

    return rates_total;
}

// n - order number of a currency pair.
void UpdateArrays(int n)
{
    ArrayFree(PairData[n].Open);
    ArrayFree(PairData[n].High);
    ArrayFree(PairData[n].Low);
    ArrayFree(PairData[n].Close);
    ArrayFree(PairData[n].Time);

    int NewBars = Bars(Pairs[n], _Period);
    // The number of requested bars (at its peak) should be the same for all pairs:
    if (NewBars < PairsLastBars[ArrayMaximum(PairsLastBars)]) NewBars = PairsLastBars[ArrayMaximum(PairsLastBars)];

    int to_copy = NewBars - PairsLastBars[n];
    if (to_copy < 2) to_copy = 2; // At least 2 bars.
    // Adjust for max amount of bars set by user.
    if ((MaxBars > 0) && (to_copy > MaxBars)) to_copy = MaxBars;
    if (to_copy > NewBars) to_copy = NewBars;
    
    int cp = CopyOpen(Pairs[n], _Period, 0, to_copy, PairData[n].Open);
    if (cp == -1) return; // No data yet.
    if (CopyHigh(Pairs[n], _Period, 0, to_copy, PairData[n].High) == -1) return;
    if (CopyLow(Pairs[n], _Period, 0, to_copy, PairData[n].Low) == -1) return;
    if (CopyClose(Pairs[n], _Period, 0, to_copy, PairData[n].Close) == -1) return;
    if (CopyTime(Pairs[n], _Period, 0, to_copy, PairData[n].Time) == -1) return;

    if (cp == to_copy) PairsLastBars[n] = NewBars; // Remember the number of processed bars for the pair only if managed to copy that number of bars.
}

//+------------------------------------------------------------------+
//| Calculates moving averages over the USDX.                        |
//+------------------------------------------------------------------+
void CalculateMAs()
{
    int i;

    // Find the first bar with values.
    for (i = 0; i < RatesTotal; i++)
        if (Open[i] != EMPTY_VALUE) break;

    if (RatesTotal - 1 - i < MA_Period1) return; // Not enough bars.
    if (MA_Period1 > 0) CalcMA(MA1, i, MA_Mode1, MA_Period1, MA_PriceType1); // Calculate only if MA period is given.
    if (RatesTotal - 1 - i < MA_Period2) return;
    if (MA_Period2 > 0) CalcMA(MA2, i, MA_Mode2, MA_Period2, MA_PriceType2);

    LastMA = RatesTotal - 1;
}

// MA calculation methods encapsulator.
void CalcMA(double &MA[], int i, ENUM_MA_METHOD ma_method, int ma_period, ENUM_APPLIED_PRICE ma_pricetype)
{
    switch (ma_method)
    {
    case MODE_SMA:
        CalcSMA(MA, i, ma_period, ma_pricetype);
        break;
    case MODE_EMA:
        CalcEMA(MA, i, ma_period, ma_pricetype);
        break;
    case MODE_SMMA:
        CalcSMMA(MA, i, ma_period, ma_pricetype);
        break;
    case MODE_LWMA:
        CalcLWMA(MA, i, ma_period, ma_pricetype);
        break;
    default:
        CalcSMA(MA, i, ma_period, ma_pricetype);
        break;
    }
}

// Simple MA Calculation.
void CalcSMA(double &MA[], int i, int ma_period, ENUM_APPLIED_PRICE ma_pricetype)
{
    // LastMA - 1 to always recalc last two bars.
    int oldest_bar_to_calc = MathMax(i + ma_period, LastMA - 1);
    // From new to old.
    for (int j = RatesTotal - 1; j >= oldest_bar_to_calc; j--)
    {
        double MA_Sum = 0;

        for (int k = j; k > j - ma_period; k--)
        {
            MA_Sum += GetPrice(ma_pricetype, k);
        }

        MA[j] = MA_Sum / ma_period;
    }
}

// Exponential MA Calculation.
void CalcEMA(double &MA[], int i, int ma_period, ENUM_APPLIED_PRICE ma_pricetype)
{
    // LastMA - 1 to always recalc last two bars.
    int oldest_bar_to_calc = MathMax(i, LastMA - 1);
    // From old to new.
    for (int j = oldest_bar_to_calc; j <= RatesTotal - 1; j++)
    {
        if (j == i) MA[j] = GetPrice(ma_pricetype, j);
        else
        {
            double coeff = 2 / (double(ma_period) + 1);
            MA[j] = GetPrice(ma_pricetype, j) * coeff + MA[j - 1] * (1 - coeff);
        }
    }
}

// Smoothed MA Calculation (Exponential with coeff = 1 / N).
void CalcSMMA(double &MA[], int i, int ma_period, ENUM_APPLIED_PRICE ma_pricetype)
{
    // LastMA - 1 to always recalc last two bars.
    int oldest_bar_to_calc = MathMax(i, LastMA - 1);
    // From old to new.
    for (int j = oldest_bar_to_calc; j <= RatesTotal - 1; j++)
    {
        if (j == i) MA[j] = GetPrice(ma_pricetype, j);
        else
        {
            double coeff = 1 / double(ma_period);
            MA[j] = GetPrice(ma_pricetype, j) * coeff + MA[j - 1] * (1 - coeff);
        }
    }
}

// Linear Weighted MA Calculation.
void CalcLWMA(double &MA[], int i, int ma_period, ENUM_APPLIED_PRICE ma_pricetype)
{
    // LastMA - 1 to always recalc last two bars.
    int oldest_bar_to_calc = MathMax(i + ma_period, LastMA - 1);
    // From new to old.
    for (int j = RatesTotal - 1; j >= oldest_bar_to_calc; j--)
    {
        double MA_Sum = 0;
        int weight = ma_period;
        int weight_sum = 0;
        for (int k = j; k > j - ma_period; k--)
        {
            weight_sum += weight;
            MA_Sum += GetPrice(ma_pricetype, k) * weight;
            weight--;
        }
        MA[j] = MA_Sum / weight_sum;
    }
}


// Get price of a given type for a given bar.
double GetPrice(ENUM_APPLIED_PRICE price_type, int n)
{
    switch(price_type)
    {
    case PRICE_OPEN:
        return Open[n];
    case PRICE_HIGH:
        return High[n];
    case PRICE_LOW:
        return Low[n];
    case PRICE_CLOSE:
        return Close[n];
    case PRICE_MEDIAN:
        return (High[n] + Low[n]) / 2;
    case PRICE_TYPICAL:
        return (High[n] + Low[n] + Close[n]) / 3;
    case PRICE_WEIGHTED:
        return (High[n] + Low[n] + Close[n] * 2) / 4;
    }

    return 0;
}
//+------------------------------------------------------------------+