//+------------------------------------------------------------------+
//|                                                         USDX.mq4 |
//|                                 Copyright © 2012-2025, EarnForex |
//|                                        https://www.earnforex.com |
//+------------------------------------------------------------------+
#property copyright "www.EarnForex.com, 2012-2025"
#property link      "https://www.earnforex.com/metatrader-indicators/USDX/"
#property version   "1.04"
#property strict

#property description "USDX - indicator of the US Dollar Index."
#property description "Displays a USDX chart in a separate window of the current chart."
#property description "Based on EUR/USD, USD/JPY, GBP/USD, USD/CAD, USD/CHF, and USD/SEK."
#property description "All these pairs should be added to Market Watch for indicator to work."
#property description "Two customizable moving averages can be applied to the index."
#property description "Can be easily modified via input parameters to calculate any currency index."

#property indicator_separate_window
#property indicator_buffers 3

#property indicator_type1  DRAW_LINE
#property indicator_color1 clrRed
#property indicator_width1 2
#property indicator_label1 "USDX"

#property indicator_type2  DRAW_LINE
#property indicator_color2 clrWhite

#property indicator_type3  DRAW_LINE
#property indicator_color3 clrLime

input ENUM_APPLIED_PRICE USDX_PriceType = PRICE_CLOSE;
input string IndexPairs = "EURUSD, USDJPY, GBPUSD, USDCAD, USDSEK, USDCHF"; // Update currency pairs' names only if different from your broker's.
input string IndexCoefficients = "-0.576, 0.136, -0.119, 0.091, 0.042, 0.036"; // IndexCoefficients: Do not change for USDX.
input double IndexInitialValue = 50.14348112; // IndexInitialValue: Do not change for USDX.
input int MA_Period1 = 13;
input int MA_Period2 = 17;
input ENUM_MA_METHOD MA_Mode1 = MODE_SMA;
input ENUM_MA_METHOD MA_Mode2 = MODE_SMA;

// Indicator buffers:
double USDX[];
double MA1[];
double MA2[];

// Global variables:
string Pairs[];
double Coefficients[];
int Pairs_NotFound[]; // 0 - OK, 1 - first failure, 2 - final failure.
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
    ArrayResize(Coefficients, count);
    ArrayResize(Pairs_NotFound, count);
    ArrayInitialize(Pairs_NotFound, 0);
    
    n = 0;
    coef_n = 0;
    int prev_n = 0, prev_coef_n = 0, return_value = 0;
    for (int i = 0; i < count; i++)
    {
        n = StringFind(IndexPairs, ",", n + 1);
        coef_n = StringFind(IndexCoefficients, ",", coef_n + 1);
        // To remove the trailing comma
        if (prev_n > 0) prev_n++;
        if (prev_coef_n > 0) prev_coef_n++;
        Pairs[i] = StringTrimRight(StringTrimLeft(StringSubstr(IndexPairs, prev_n, n - prev_n)));
        Coefficients[i] = StrToDouble(StringSubstr(IndexCoefficients, prev_coef_n, coef_n - prev_coef_n));
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

    return return_value;
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
    
    SetIndexBuffer(0, USDX);
    SetIndexEmptyValue(0, EMPTY_VALUE);

    if (MA_Period1 > 0)
    {
        SetIndexBuffer(1, MA1);
        SetIndexLabel(1, "MA(" + IntegerToString(MA_Period1) + ")");
        SetIndexEmptyValue(1, 0);
    } else SetIndexStyle(1, DRAW_NONE);

    if (MA_Period2 > 0)
    {
        SetIndexBuffer(2, MA2);
        SetIndexLabel(2, "MA(" + IntegerToString(MA_Period2) + ")");
        SetIndexEmptyValue(2, 0);
    } else SetIndexStyle(2, DRAW_NONE);

    IndicatorShortName("USDX");

    return INIT_SUCCEEDED;
}

int limit;

int OnCalculate(const int rates_total,
                const int prev_calculated,
                const datetime &time[],
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
    if (!AllFound) return 0;

    int counted_bars = prev_calculated;
    if (counted_bars > 0) counted_bars--;
    limit = Bars - counted_bars;
    for (int i = 0; i < limit; i++)
    {
        USDX[i] = IndexInitialValue;
        for (int j = 0; j < ArraySize(Pairs); j++)
        {
            int index = iBarShift(Pairs[j], 0, Time[i], false);
            if (index < 0) continue;
            if (iBars(Pairs[j], 0) == 0) return 0; // Chart history not loaded yet.
            switch (USDX_PriceType)
            {
            case PRICE_OPEN:
                USDX[i] *= MathPow(iOpen(Pairs[j], 0, index), Coefficients[j]);
                break;
            case PRICE_HIGH:
                USDX[i] *= MathPow(iHigh(Pairs[j], 0, index), Coefficients[j]);
                break;
            case PRICE_LOW:
                USDX[i] *= MathPow(iLow(Pairs[j], 0, index), Coefficients[j]);
                break;
            case PRICE_CLOSE:
                USDX[i] *= MathPow(iClose(Pairs[j], 0, index), Coefficients[j]);
                break;
            case PRICE_MEDIAN:
                USDX[i] *= MathPow((iHigh(Pairs[j], 0, index) + iLow(Pairs[j], 0, index)) / 2, Coefficients[j]);
                break;
            case PRICE_TYPICAL:
                USDX[i] *= MathPow((iHigh(Pairs[j], 0, index) + iLow(Pairs[j], 0, index) + iClose(Pairs[j], 0, index)) / 3, Coefficients[j]);
                break;
            case PRICE_WEIGHTED:
                USDX[i] *= MathPow((iHigh(Pairs[j], 0, index) + iLow(Pairs[j], 0, index) + iClose(Pairs[j], 0, index) * 2) / 4, Coefficients[j]);
                break;
            default:
                USDX[i] *= MathPow(iClose(Pairs[j], 0, index), Coefficients[j]);
                break;
            }
        }
    }

    CalculateMAs();

    return rates_total;
}

//+------------------------------------------------------------------+
//| Calculates moving averages over the USDX.                        |
//+------------------------------------------------------------------+
void CalculateMAs()
{
    if (Bars < MA_Period1) return; // Not enough bars.
    if (MA_Period1 > 0) CalcMA(MA1, MA_Mode1, MA_Period1); // Calculate only if an MA period is given.
    if (Bars < MA_Period2) return;
    if (MA_Period2 > 0) CalcMA(MA2, MA_Mode2, MA_Period2);
}

// MA calculation methods encapsulator.
void CalcMA(double &MA[], int ma_method, int ma_period)
{
    switch (ma_method)
    {
    case MODE_SMA:
        CalcSMA(MA, ma_period);
        break;
    case MODE_EMA:
        CalcEMA(MA, ma_period);
        break;
    case MODE_SMMA:
        CalcSMMA(MA, ma_period);
        break;
    case MODE_LWMA:
        CalcLWMA(MA, ma_period);
        break;
    default:
        CalcSMA(MA, ma_period);
        break;
    }
}

// Simple MA Calculation.
void CalcSMA(double &MA[], int ma_period)
{
    int local_limit = limit;
    if (local_limit > Bars - ma_period) local_limit = Bars - ma_period;
    // From new to old
    for (int j = 0; j < local_limit; j++)
    {
        double MA_Sum = 0;

        for (int k = j; k < j + ma_period; k++)
        {
            MA_Sum += USDX[k];
        }

        MA[j] = MA_Sum / ma_period;
    }
}

// Exponential MA Calculation.
// ma_period is double here because int produces wrong division results in coeff calculation below.
void CalcEMA(double &MA[], double ma_period)
{
    // From old to new.
    for (int j = limit - 1; j >= 0; j--)
    {
        // Cannot calculate EMA for empty value.
        if (USDX[j] == EMPTY_VALUE) continue;

        // If no previous MA value, take the price value.
        if ((MA[j + 1] == 0) || (MA[j + 1] == EMPTY_VALUE)) MA[j] = USDX[j];
        else
        {
            double coeff = 2 / (ma_period + 1);
            MA[j] = USDX[j] * coeff + MA[j + 1] * (1 - coeff);
        }
    }
}

// Smoothed MA Calculation (Exponential with coeff = 1 / N).
// ma_period is double here because int produces wrong division results in coeff calculation below.
void CalcSMMA(double &MA[], double ma_period)
{
    // From old to new.
    for (int j = limit - 1; j >= 0; j--)
    {
        // Cannot calculate SMMA for empty value.
        if (USDX[j] == EMPTY_VALUE) continue;

        // If no previous MA value, take the price value.
        if ((MA[j + 1] == 0) || (MA[j + 1] == EMPTY_VALUE)) MA[j] = USDX[j];
        else
        {
            double coeff = 1 / ma_period;
            MA[j] = USDX[j] * coeff + MA[j + 1] * (1 - coeff);
        }
    }
}

// Linear Weighted MA Calculation.
void CalcLWMA(double &MA[], int ma_period)
{
    int local_limit = limit;
    if (local_limit > Bars - ma_period) local_limit = Bars - ma_period;
    // From new to old.
    for (int j = 0; j < local_limit; j++)
    {
        double MA_Sum = 0;
        int weight = ma_period;
        int weight_sum = 0;
        for (int k = j; k < j + ma_period; k++)
        {
            weight_sum += weight;
            MA_Sum += USDX[k] * weight;
            weight--;
        }
        MA[j] = MA_Sum / weight_sum;
    }
}
//+------------------------------------------------------------------+