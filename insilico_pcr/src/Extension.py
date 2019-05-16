# Author: Izaak Coleman
# email: izaak.coleman1@gmail.com

__version__ = 0.1

# C++-like enums to control extension direction
class Extend(object): pass
class Forward(Extend): pass
class Backward(Extend): pass

class Extension:
  """Encodes a sequence present in one or more fastq databases.
     
     The encoded sequence (extension) is constructed from the iPCR module which manages
     construction of the extension. The class also encodes which sequencing database the
     extension is found in, and for each database the invariant U_q(p1) - U_m(p1) that 
     must be maintained throughout construction of the extension."""

  def __init__(self, probe):
    """probe: An extension of an Extension object begins with an intial probe."""
    self.extension = probe # initialise the extension
    self.databases = dict() # Will store a dictionary of k = database_id, v = U_q(p1) - U_m(p1) invariant for database
    self.init_probe_len = len(probe) # required to cut probe from self.extension for backward extension

  def extend(self, base, direction):
    """self.extension will be extended in the forward direction (towards p2) as request by iPCR."""
    if base not in 'ATCG':
      raise Exception(f"input base {base} is not a valid base to extend Extension, base must be [ATCG]")
    if direction == Forward:
      self.extension = self.extension + base
    elif direction == Backward:
      self.extension = base + self.extension
    else:
      raise Exception(f"Direction {direction} is an invalid direction for extension")

  def cut_initializing_probe(self):
    """ Cleaves the initializing probe from self.extension.

        Once extension has reached p1 + sigma^k + p2*, it is now desired to remove the intial probe, 
        and extend in the backward direction from sigma^k + p2* to recover the final extension
        p1* + sigma^k + p2*. In order to do this , the inital input probe must be cleaved."""
    self.extension = self.extension[self.init_probe_len:]

  def __len__(self):
    """Consider len(Extension) to be len(self.extension)."""
    return len(self.extension)

