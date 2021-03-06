// Copyright 2017 GoogleGenomics R package authors. All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the 'License');
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an 'AS IS' BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <memory>
#include <string>

#include <grpc++/channel.h>
#include <grpc++/client_context.h>
#include <grpc++/support/status.h>
#include <grpc++/support/sync_stream.h>
#include <grpc/grpc.h>
#include <google/genomics/v1/variants.grpc.pb.h>
#include <google/genomics/v1/variants.pb.h>

#include <client-common.h>

using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

using ::grpc::ClientContext;
using ::grpc::ClientReader;
using ::grpc::Status;

using namespace ::google::genomics::v1;
using namespace ::googlegenomics;

extern "C" {
STREAMED_METHOD(StreamingVariantService, StreamVariants, StreamVariantsRequest,
                StreamVariantsResponse)

METHOD(VariantServiceV1, ImportVariants, ImportVariantsRequest,
       google::longrunning::Operation)

METHOD(VariantServiceV1, CreateVariantSet, CreateVariantSetRequest, VariantSet)

METHOD(VariantServiceV1, ExportVariantSet, ExportVariantSetRequest,
       google::longrunning::Operation)

METHOD(VariantServiceV1, GetVariantSet, GetVariantSetRequest, VariantSet)

METHOD(VariantServiceV1, SearchVariantSets, SearchVariantSetsRequest,
       SearchVariantSetsResponse)

METHOD(VariantServiceV1, DeleteVariantSet, DeleteVariantSetRequest,
       google::protobuf::Empty)

METHOD(VariantServiceV1, UpdateVariantSet, UpdateVariantSetRequest, VariantSet)

METHOD(VariantServiceV1, SearchVariants, SearchVariantsRequest,
       SearchVariantsResponse)

METHOD(VariantServiceV1, CreateVariant, CreateVariantRequest, Variant)

METHOD(VariantServiceV1, UpdateVariant, UpdateVariantRequest, Variant)

METHOD(VariantServiceV1, DeleteVariant, DeleteVariantRequest,
       google::protobuf::Empty)

METHOD(VariantServiceV1, GetVariant, GetVariantRequest, Variant)

METHOD(VariantServiceV1, MergeVariants, MergeVariantsRequest,
       google::protobuf::Empty)

METHOD(VariantServiceV1, SearchCallSets, SearchCallSetsRequest,
       SearchCallSetsResponse)

METHOD(VariantServiceV1, CreateCallSet, CreateCallSetRequest, CallSet)

METHOD(VariantServiceV1, UpdateCallSet, UpdateCallSetRequest, CallSet)

METHOD(VariantServiceV1, DeleteCallSet, DeleteCallSetRequest,
       google::protobuf::Empty)

METHOD(VariantServiceV1, GetCallSet, GetCallSetRequest, CallSet)
}  // extern C
