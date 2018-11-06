/*
Module : AAChineseCalendar.cpp
Purpose: Implementation for the algorithms which convert between the Gregorian and Julian calendars and the Chinese calendar
Created: Rui Xia / 04-11-2018

Copyright (c) 2004 - 2018 by PJ Naughter (Web: www.naughter.com, Email: pjna@naughter.com)

All rights reserved.

Copyright / Usage Details:

You are allowed to include the source code in any product (commercial, shareware, freeware or otherwise) 
when your product is released in binary form. You are allowed to modify the source code in any way you want 
except you cannot modify the copyright details at the top of each module. If you want to distribute source 
code with your application, then you are only allowed to distribute versions released by the author. This is 
to maintain a single distribution point for the source code. 

*/


//////////////////////////// Includes /////////////////////////////////////////

#include "stdafx.h"
#include "AAChineseCalendar.h"
#include "AACoordinateTransformation.h"
#include "AASun.h"
#include "AAMoon.h"
#include "AAEarth.h"
#include "AANutation.h"
#include "AAMoonIlluminatedFraction.h"
#include "AAELPMPP02.h"
#include "AAELP2000.h"
#include "AAMoonPhases.h"
#include "AAPrecession.h"
#include "AAAberration.h"
#include "AAFK5.h"
#include "AADynamicalTime.h"
#include <cmath>
#include <unordered_map>
using namespace std;


//////////////////////////// Implementation ///////////////////////////////////

const std::unordered_set<CAAChineseCalendar::SolarTerm> CAAChineseCalendar::AllSolarTerms{
  SolarTerm::SpringCommences,
  SolarTerm::SpringShowers,
  SolarTerm::InsectsWaken,
  SolarTerm::VernalEquinox,
  SolarTerm::BrightAndClear,
  SolarTerm::CornRain,
  SolarTerm::SummerCommences,
  SolarTerm::CornForms,
  SolarTerm::CornOnEar,
  SolarTerm::SummerSolstice,
  SolarTerm::ModerateHeat,
  SolarTerm::GreatHeat,
  SolarTerm::AutumnCommences,
  SolarTerm::EndOfHeat,
  SolarTerm::WhiteDew,
  SolarTerm::AutumnalEquinox,
  SolarTerm::ColdDew,
  SolarTerm::Frost,
  SolarTerm::WinterCommences,
  SolarTerm::LightSnow,
  SolarTerm::HeavySnow,
  SolarTerm::WinterSolstice,
  SolarTerm::ModerateCold,
  SolarTerm::SevereCold,
};

const std::unordered_set<CAAChineseCalendar::SolarTerm> CAAChineseCalendar::MidpointSolarTerms{
  SolarTerm::SpringShowers,
  SolarTerm::VernalEquinox,
  SolarTerm::CornRain,
  SolarTerm::CornForms,
  SolarTerm::SummerSolstice,
  SolarTerm::GreatHeat,
  SolarTerm::EndOfHeat,
  SolarTerm::AutumnalEquinox,
  SolarTerm::Frost,
  SolarTerm::LightSnow,
  SolarTerm::WinterSolstice,
  SolarTerm::SevereCold,
};

double CAAChineseCalendar::MoonElongation(double JD) {
  double sunDistance = 149597870.7 * CAAEarth::RadiusVector(JD, true);
  double moonDistance = CAAMoon::RadiusVector(JD);

  CAA3DCoordinate moonDerivative;
  CAA3DCoordinate moon = CAAELPMPP02::EclipticRectangularCoordinatesJ2000(JD, CAAELPMPP02::Correction::LLR, &moonDerivative);
  double moonEclipticLongitude = CAACoordinateTransformation::RadiansToDegrees(std::atan2(moon.Y, moon.X));
  double moonEclipticLatitude = CAACoordinateTransformation::RadiansToDegrees(std::asin(moon.Z / std::sqrt(moon.X * moon.X + moon.Y * moon.Y + moon.Z * moon.Z)));

  // Precess coordinate from J2000 to the dynamical ecliptic and equinox of the date defined in VSOP theory
  CAA2DCoordinate moon2d = CAAPrecession::PrecessEcliptic(moonEclipticLongitude, moonEclipticLatitude, 2451545.0, JD);
  moonEclipticLongitude = moon2d.X;
  moonEclipticLatitude = moon2d.Y;

  // VSOP to FK5
  moonEclipticLongitude += CAAFK5::CorrectionInLongitude(moonEclipticLongitude, moonEclipticLatitude, JD);
  moonEclipticLatitude += CAAFK5::CorrectionInLatitude(moonEclipticLongitude, JD);

  // Nutation correction
  moonEclipticLongitude += CAACoordinateTransformation::DMSToDegrees(0, 0, CAANutation::NutationInLongitude(JD));

  // Aberration correction
  moonEclipticLongitude -= (0.005775518 * CAAMoon::RadiusVector(JD) / 149597870.7 * (360 / 27.3));

  // Apparent coordinate of sun (FK5 + nutation + abberation)
  double sunEclipticLongitude = CAASun::ApparentEclipticLongitude(JD, true);

  return CAACoordinateTransformation::MapTo0To360Range(moonEclipticLongitude - sunEclipticLongitude + 180) - 180;
}

long CAAChineseCalendar::DayFromJulian(double JD, double timezone) {
  JD = CAADynamicalTime::TT2UTC(JD) + timezone / 24;
  return CAADate::INT(JD + 0.5 + 1e-8);
}

double CAAChineseCalendar::JulianFromDay(long day, double timezone) {
  return CAADynamicalTime::UTC2TT((double)day - 0.5 - timezone / 24 + 1e-8);
}

double CAAChineseCalendar::PreviousNewMoon(double JD, bool inclusive) {
  double elongation = MoonElongation(JD);
  JD -= std::fmod(elongation + 360 + (inclusive ? 1 : -1) * 360 / (27.3 * 24 * 3600), 360) / 360 * 27.3;

  double Correction = 0;
  do {
    elongation = MoonElongation(JD);
    Correction = -elongation / 360 * 27.3;
    JD += Correction;
  } while (fabs(Correction) > 0.000001); //Corresponds to an error of 0.0864 of a second
  return JD;
}

double CAAChineseCalendar::NextNewMoon(double JD, bool inclusive) {
  double elongation = MoonElongation(JD);
  JD += std::fmod(-elongation + 360 + (inclusive ? 1 : -1) * 360 / (27.3 * 24 * 3600), 360) / 360 * 27.3;

  double Correction = 0;
  do {
    elongation = MoonElongation(JD);
    Correction = -elongation / 360 * 27.3;
    JD += Correction;
  } while (fabs(Correction) > 0.000001); //Corresponds to an error of 0.0864 of a second
  return JD;
}

long CAAChineseCalendar::PreviousNewMoon(long day, double timezone, bool inclusive) {
  double JD = JulianFromDay(day + (inclusive ? 1 : 0), timezone);
  JD = PreviousNewMoon(JD);
  return DayFromJulian(JD, timezone);
}

long CAAChineseCalendar::NextNewMoon(long day, double timezone, bool inclusive) {
  double JD = JulianFromDay(day + (inclusive ? 0 : 1), timezone);
  JD = NextNewMoon(JD);
  return DayFromJulian(JD, timezone);
}

double CAAChineseCalendar::PreviousSolarTerm(double JD, const std::unordered_set<SolarTerm> &solarTerms, bool inclusive) {
  double longitude = CAASun::ApparentEclipticLongitude(JD, true);
  std::vector<std::pair<double, double>> offsets;
  double delta = (inclusive ? 1 : -1) * 360 / (365.24219 * 24 * 3600);
  for (SolarTerm solarTerm: solarTerms) {
    int indexFromNorthwardEquinox = static_cast<int>(solarTerm) - static_cast<int>(SolarTerm::VernalEquinox);
    double solarTermLongitude = 15 * indexFromNorthwardEquinox;
    double offset = CAACoordinateTransformation::MapTo0To360Range(longitude - solarTermLongitude + delta);
    offsets.push_back(std::make_pair(offset, solarTermLongitude));
  }
  auto minOffsetIterator = std::min_element(offsets.begin(), offsets.end(),
      [](const std::pair<double, double> &a, const std::pair<double, double> &b) {
    return a.first < b.first;
  });
  JD -= minOffsetIterator->first / 360 * 365.24219;

  double Correction = 0;
  do
  {
    const double longitude = CAASun::ApparentEclipticLongitude(JD, true);
    Correction = 58 * sin(CAACoordinateTransformation::DegreesToRadians(minOffsetIterator->second - longitude)); 
    JD += Correction;
  }
  while (fabs(Correction) > 0.000001); //Corresponds to an error of 0.0864 of a second

  return JD;
}

double CAAChineseCalendar::NextSolarTerm(double JD, const std::unordered_set<SolarTerm> &solarTerms, bool inclusive) {
  double longitude = CAASun::ApparentEclipticLongitude(JD, true);
  std::vector<std::pair<double, double>> offsets;
  double delta = (inclusive ? 1 : -1) * 360 / (365.24219 * 24 * 3600);
  for (SolarTerm solarTerm: solarTerms) {
    int indexFromNorthwardEquinox = static_cast<int>(solarTerm) - static_cast<int>(SolarTerm::VernalEquinox);
    double solarTermLongitude = 15 * indexFromNorthwardEquinox;
    double offset = CAACoordinateTransformation::MapTo0To360Range(solarTermLongitude - longitude + delta);
    offsets.push_back(std::make_pair(offset, solarTermLongitude));
  }
  auto minOffsetIterator = std::min_element(offsets.begin(), offsets.end(),
      [](const std::pair<double, double> &a, const std::pair<double, double> &b) {
    return a.first < b.first;
  });
  JD += minOffsetIterator->first / 360 * 365.24219;

  double Correction = 0;
  do
  {
    const double longitude = CAASun::ApparentEclipticLongitude(JD, true);
    Correction = 58 * sin(CAACoordinateTransformation::DegreesToRadians(minOffsetIterator->second - longitude)); 
    JD += Correction;
  }
  while (fabs(Correction) > 0.000001); //Corresponds to an error of 0.0864 of a second

  return JD;
}

std::vector<long> CAAChineseCalendar::NewMoonDays(long start, long end, double timezone) {
  std::vector<long> result;
  double endJD = JulianFromDay(end, timezone) + 0.01;
  for (long day = NextNewMoon(start, timezone, false);
      JulianFromDay(day, timezone) <= endJD;
      day = NextNewMoon(day, timezone, false)) {
    result.push_back(day);
  }
  return result;
}

int CAAChineseCalendar::CalculateLeapMonth(const std::vector<long> &newMoonDays, double timezone) {
  if (newMoonDays.size() == 12) {
    return -1;
  }

  int lunarMonth = 12;
  double termJD = NextSolarTerm(JulianFromDay(newMoonDays[0], timezone), MidpointSolarTerms, true);
  for (int i = 0; i < 12; ++i) { // 12 mid-term climates
    if (termJD < JulianFromDay(newMoonDays[i + 1], timezone)) {
      termJD = NextSolarTerm(termJD, MidpointSolarTerms, false);
      lunarMonth = lunarMonth % 12 + 1;
    } else {
      // Found month without mid-term climate
      return (lunarMonth + 10) % 12 + 1;
    }
  }
  return (lunarMonth + 10) % 12 + 1;
}

CAAChineseCalendar::CalendarDate CAAChineseCalendar::JulianToChinese(long Year, long Month, long Day, bool bGregorianCalendar) noexcept {
  double timezone = +8;
  CAADate date(Year, Month, Day, bGregorianCalendar);
  double JD = JulianFromDay(CAADate::INT(CAADynamicalTime::TT2UTC(date.Julian() + 0.5)), timezone);
  double winterSolsticeJD = PreviousSolarTerm(JD + 1, {SolarTerm::WinterSolstice}, false);
  double nextWinterSolsticeJD = NextSolarTerm(JD + 1, {SolarTerm::WinterSolstice}, true);
  double summerSolsticeJD = NextSolarTerm(CAADate(Year, 1, 1, true).Julian(), {SolarTerm::SummerSolstice});
  int yearOffset = winterSolsticeJD > summerSolsticeJD ? 0 : -1;

  // Get new moon dates
  static std::unordered_map<long, std::pair<long, std::vector<long>>> newMoonBuffer;
  long startDay = DayFromJulian(winterSolsticeJD, timezone);
  long endDay = DayFromJulian(nextWinterSolsticeJD, timezone);
  std::vector<long> newMoons;
  long monthStart;

  if (newMoonBuffer.count(startDay)) {
    monthStart = newMoonBuffer[startDay].first;
    newMoons = newMoonBuffer[startDay].second;
  } else {
    newMoons = CAAChineseCalendar::NewMoonDays(startDay, endDay, timezone);
    monthStart = PreviousNewMoon(newMoons[0], timezone, false);
    newMoonBuffer[startDay] = std::make_pair(monthStart, newMoons);
  }
  int leapMonth = CAAChineseCalendar::CalculateLeapMonth(newMoons, timezone);

  int month = 11;
  bool leaped = false;
  for (int i = 0; JD + 0.01 >= JulianFromDay(newMoons[i], timezone) && i < newMoons.size(); ++i) {
    if (!leaped && month == leapMonth) {
      leaped = true;
    } else {
      leaped = false;
      month = month % 12 + 1;
      if (month == 1) {
        yearOffset++;
      }
    }
    monthStart = newMoons[i];
  }
  CalendarDate lunarDate{Year + yearOffset, month, std::lround(JD - JulianFromDay(monthStart, timezone)) + 1, leaped && month == leapMonth};
  return lunarDate;
}