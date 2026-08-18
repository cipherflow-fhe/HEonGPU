// Copyright 2025-2026 Yanbin Li
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
// Developer: Yanbin Li

#ifndef HEONGPU_METADATA_H
#define HEONGPU_METADATA_H

namespace heongpu
{
    struct Metadata // @company CipherFlow
    {
        bool is_ringt;
        int level;
        double scale;
        int log_slot_count;
    };

} // namespace heongpu
#endif // HEONGPU_METADATA_H
