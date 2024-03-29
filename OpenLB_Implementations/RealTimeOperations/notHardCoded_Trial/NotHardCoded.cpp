/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2006-2019 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

#include "NotHardCoded.h"
#include "NotHardCoded.hh"

using T = FLOATING_POINT_TYPE;

std::vector<std::string> stringStart;
std::vector<std::shared_ptr<IndicatorF3D<T>>> indicators;

int main( int argc, char* argv[] )
{


  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  clout << "Not Hard Coded Model..." << std::endl;

  redisContext *redis = redisConnect("127.0.0.1", 6379);

  redisReply *response = (redisReply *)redisCommand(redis, "FLUSHALL");
  freeReplyObject(response);
  redisFree(redis);

  bool RedisStartValues = true;

  while(RedisStartValues)
  {
    redisContext *redis = redisConnect("127.0.0.1", 6379);
    redisReply *reply = (redisReply *)redisCommand(redis, "LRANGE converterValues 0 -1");

    if (reply != NULL && reply->type == REDIS_REPLY_ARRAY && reply->elements > 7) 
    {
      for (size_t i = 0; i < reply->elements; ++i) {
          stringStart.push_back(reply->element[i]->str);
      }

      if(stringStart.size() > 6)
      {
        RedisStartValues = false;

        std::cout << "Start Properties: " << stringStart << std::endl;
      }
      freeReplyObject(reply);
    } 
    else 
    {
        std::cout << "Failed to retrieve strings from Redis." << std::endl;
    }

    redisFree(redis);
    std::this_thread::sleep_for(std::chrono::seconds(2));
  }

  std::cout << "Creating Model..." << std::endl;
  ModelMSO<T> model(stringStart);
  std::cout << "Model Created..." << std::endl;
  
  bool getGeometryAndBoundaryData = false;
  while(!getGeometryAndBoundaryData)
  {
    getGeometryAndBoundaryData = model.getGeometryData();
    std::this_thread::sleep_for(std::chrono::seconds(2));
  }

  getGeometryAndBoundaryData = false;
  while(!getGeometryAndBoundaryData)
  {
    getGeometryAndBoundaryData = model.getAndSetBoundaryData();
    std::this_thread::sleep_for(std::chrono::seconds(2));
  }

  std::cout << "Setting Geometry Data..." << std::endl;
  model.setGeometryData();

  std::cout << "Starting Simulation..." << std::endl;
  model.simulate();
}
