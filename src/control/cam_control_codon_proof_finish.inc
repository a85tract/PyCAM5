       if (cam_control_tag_out /= int(CAM_CONTROL_PROOF_TAG, c_int64_t)) then
          write(iulog,*) CAM_CONTROL_PROOF_LABEL//' Codon tag roundtrip failed'
          stop 2
       endif
       write(iulog,*) CAM_CONTROL_PROOF_LABEL//' implementation = codon'
       cam_control_proof_seen = .true.
    endif
