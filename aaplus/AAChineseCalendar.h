/*
Module : AAMoslemCalendar.h
Purpose: Implementation for the algorithms which convert between the Julian and Moslem calendars
Created: PJN / 04-02-2004

Copyright (c) 2004 - 2018 by PJ Naughter (Web: www.naughter.com, Email: pjna@naughter.com)

All rights reserved.

Copyright / Usage Details:

You are allowed to include the source code in any product (commercial, shareware, freeware or otherwise)
when your product is released in binary form. You are allowed to modify the source code in any way you want
except you cannot modify the copyright details at the top of each module. If you want to distribute source
code with your application, then you are only allowed to distribute versions released by the author. This is
to maintain a single distribution point for the source code.

*/


/////////////////////// Macros / Defines //////////////////////////////////////

#if _MSC_VER > 1000
#pragma once
#endif //#if _MSC_VER > 1000

#ifndef __AACHINESECALENDAR_H__
#define __AACHINESECALENDAR_H__

#ifndef AAPLUS_EXT_CLASS
#define AAPLUS_EXT_CLASS
#endif //#ifndef AAPLUS_EXT_CLASS


/////////////////////// Includes //////////////////////////////////////////////

#include "AADate.h"
#include <vector>
#include <algorithm>
#include <unordered_set>

/////////////////////// Classes ///////////////////////////////////////////////

class AAPLUS_EXT_CLASS CAAChineseCalendar
{
public:
  // Solar terms: https://en.wikipedia.org/wiki/Solar_term
  // English names taken from: http://www.weather.gov.hk/gts/time/24solarterms.htm
  // Longtitude of the sun can be calculated as (int)SolarTerm * 15 in degrees.
  enum class SolarTerm {
    SpringCommences = 0,    // Lichun
    SpringShowers = 1,      // Yushui
    InsectsWaken = 2,       // Jingzhe
    VernalEquinox = 3,      // Chunfen
    BrightAndClear = 4,     // Qingming
    CornRain = 5,           // Guyu
    SummerCommences = 6,    // Lixia
    CornForms = 7,          // Xiaoman
    CornOnEar = 8,          // Mangzhong
    SummerSolstice = 9,     // Xiaozhi
    ModerateHeat = 10,      // Xiaoshu
    GreatHeat = 11,         // Dashu
    AutumnCommences = 12,   // Liqiu
    EndOfHeat = 13,         // Chushu
    WhiteDew = 14,          // Bailu
    AutumnalEquinox = 15,   // Qiufen
    ColdDew = 16,           // Hanlu
    Frost = 17,             // Shuangjiang
    WinterCommences = 18,   // Lidong
    LightSnow = 19,         // Xiaoxue
    HeavySnow = 20,         // Daxue
    WinterSolstice = 21,    // Dongzhi
    ModerateCold = 22,      // Xiaohan
    SevereCold = 23,        // Dahan
  };

  static const std::unordered_set<SolarTerm> AllSolarTerms;
  static const std::unordered_set<SolarTerm> MidpointSolarTerms;

  // Heavenly Stems (Tiangan)
  // https://en.wikipedia.org/wiki/Heavenly_Stems
  enum class HeavenlyStem {
    Jia = 0,
    Yi = 1,
    Bing = 2,
    Ding = 3,
    Wu = 4,
    Ji = 5,
    Geng = 6,
    Xin = 7,
    Ren = 8,
    Gui = 9,
  };

  // Earthly Branches (Dizhi)
  // https://en.wikipedia.org/wiki/Earthly_Branches
  enum class EarthlyBranch {
    Zi = 0,
    Chou = 1,
    Yin = 2,
    Mao = 3,
    Chen = 4,
    Si = 5,
    Wu = 6,
    Wei = 7,
    Shen = 8,
    You = 9,
    Xu = 10,
    Hai = 11,
  };

  class SexagenaryCycle {
  public:
    SexagenaryCycle();
    static std::pair<HeavenlyStem, EarthlyBranch> FromOrdinal(int ordinal);
    int ToOrdinal(HeavenlyStem  heanvenlyStem, EarthlyBranch earthlyBranch);
  };

  enum class MonthPolicy {
    Empirical,    // Pingshuo
    Astronomical, // Dingshuo -- this is the one used since Shoushi Li (Yuan)
  };

  enum class SolarTermPolicy {
    Empirical,    // Pingqi
    Astronomical, // Dingqi -- this is the one used today
  };

  struct MonthPolicyDetails {
    double daysInYear;
    double daysInMonth;

    // Period of Zhangdong, called Doufencha
    double daysInNutation;
    // A negative value indicating the length difference of tropical year over time. Applied since Tongtian Li (Song)
    double daysInYearDifference;
  };

private:
  // Get the apparent ecliptic longitude difference between the Moon and the Sun.
  static double MoonElongation(double JD);

  // Convert Julian Date to a integer representing number of days since JD0 in UTC timezone.
  static long DayFromJulian(double JD, double timezone);

  // Convert the number of days since JD0 in UTC timezone to Julian Date
  static double JulianFromDay(long day, double timezone);

  // Return the new moon date prior to the given date
  static double PreviousNewMoon(double JD, bool inclusive = false);

  // Return the new moon date after the given date
  static double NextNewMoon(double JD, bool inclusive = false);

  // Return the new moon calendar date prior to the given calendar date
  static long PreviousNewMoon(long day, double timezone, bool inclusive = false);

  // Return the new moon date after the given date
  static long NextNewMoon(long day, double timezone, bool inclusive = false);

  // Return the date of the specific Solar Term prior to the given date
  static double PreviousSolarTerm(double JD, const std::unordered_set<SolarTerm> &solarTerms = AllSolarTerms, bool inclusive = false);

  // Return the date of the specific Solar Term after the given date
  static double NextSolarTerm(double JD, const std::unordered_set<SolarTerm> &solarTerms = AllSolarTerms, bool inclusive = false);

  // Get all new moon dates between start (exclude) and end (include)
  // Used in leap months calculation
  static std::vector<long> NewMoonDays(long start, long end, double timezone);

  // Use the return value of NewMoonDates to calculate the leap month.
  // If there is 12 dates in newMoonDates, there is no leap month. The function returns -1 in the case.
  // If there is 13 dates in newMonnDaets, there one leap month after the first month which does not have a Mid-climate Term.
  static int CalculateLeapMonth(const std::vector<long> &newMoonDays, double timezone);


public:
  struct CalendarDate {
    int Year;
    int Month;
    int Day;
    bool LeapMonth;
  };

  static CalendarDate JulianToChinese(long Year, long Month, long Day, bool bGregorianCalendar) noexcept;
  // static CAACalendarDate ChineseToJulian(long Year, long Month, long Day, bool LeapMonth) noexcept;
};


#endif //#ifndef __AACHINESECALENDAR_H__
