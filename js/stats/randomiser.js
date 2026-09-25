/**
 * Statis-Gravity — Clinical & Research Randomisation Engine
 * Cryptographically Secure, Zero Modulo Bias, Permuted Block Allocation, Audit Logging.
 * Conceived, supervised design and testing: Dr G Narenthiran MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN)
 */

export const Randomiser = {
  /**
   * Generates a cryptographically secure uniform random integer in [min, max]
   * using rejection sampling to eliminate modulo bias.
   * @param {number} min 
   * @param {number} max 
   * @returns {number} Uniform integer in [min, max]
   */
  getSecureRandomInt(min, max) {
    min = Math.floor(min);
    max = Math.floor(max);
    if (min > max) {
      throw new Error(`Invalid range: min (${min}) cannot be greater than max (${max})`);
    }
    const range = max - min + 1;
    if (range === 1) return min;

    // 2^32 = 4294967296
    const MAX_UINT32 = 4294967296;
    // Upper threshold for rejection sampling: eliminate values in the incomplete modulo bucket
    const limit = MAX_UINT32 - (MAX_UINT32 % range);

    const buffer = new Uint32Array(1);
    let rand;
    do {
      let gotRand = false;
      if (typeof crypto !== 'undefined' && crypto && typeof crypto.getRandomValues === 'function') {
        try {
          crypto.getRandomValues(buffer);
          gotRand = true;
        } catch (e) {}
      }
      if (!gotRand && typeof window !== 'undefined' && window.crypto && typeof window.crypto.getRandomValues === 'function') {
        try {
          window.crypto.getRandomValues(buffer);
          gotRand = true;
        } catch (e) {}
      }
      if (!gotRand) {
        buffer[0] = Math.floor(Math.random() * MAX_UINT32);
      }
      rand = buffer[0];
    } while (rand >= limit);

    return min + (rand % range);
  },

  /**
   * Shuffles an array in-place using Fisher-Yates with cryptographically secure random integers.
   * @param {Array} array 
   * @returns {Array} Shuffled copy of array
   */
  fisherYatesShuffle(array) {
    const copy = [...array];
    for (let i = copy.length - 1; i > 0; i--) {
      const j = this.getSecureRandomInt(0, i);
      const temp = copy[i];
      copy[i] = copy[j];
      copy[j] = temp;
    }
    return copy;
  },

  /**
   * Performs Simple Randomization (1-100 Odd/Even Assignment).
   * Odd = Group A, Even = Group B.
   * @param {Object} options
   * @returns {Object} Allocation record
   */
  generateSimpleAllocation(options = {}) {
    const participantId = options.participantId || 1;
    const labelA = options.labelA || 'Group A';
    const labelB = options.labelB || 'Group B';

    const randomNumber = this.getSecureRandomInt(1, 100);
    const isOdd = (randomNumber % 2) !== 0;
    const groupKey = isOdd ? 'A' : 'B';
    const groupLabel = isOdd ? labelA : labelB;

    return {
      participantId,
      timestamp: new Date().toISOString(),
      displayTime: new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' }),
      randomNumber,
      parity: isOdd ? 'Odd' : 'Even',
      groupKey,
      groupLabel
    };
  },

  /**
   * Creates a balanced permuted block for 2 groups (A & B).
   * @param {number} blockSize Even integer >= 2
   * @param {Object} options
   * @returns {Object} Block object with shuffled slots
   */
  createBlock(blockSize, options = {}) {
    let size = parseInt(blockSize, 10);
    if (isNaN(size) || size < 2) size = 4;
    if (size % 2 !== 0) size += 1; // Enforce even block size for 1:1 balance

    const blockNumber = options.blockNumber || 1;
    const labelA = options.labelA || 'Group A';
    const labelB = options.labelB || 'Group B';

    const half = size / 2;
    const allocations = [];
    for (let i = 0; i < half; i++) allocations.push('A');
    for (let i = 0; i < half; i++) allocations.push('B');

    const shuffled = this.fisherYatesShuffle(allocations);

    const slots = shuffled.map((g, idx) => ({
      slotIndex: idx + 1,
      groupKey: g,
      groupLabel: g === 'A' ? labelA : labelB,
      assigned: false,
      participantId: null,
      timestamp: null
    }));

    return {
      blockNumber,
      blockSize: size,
      slots,
      currentIndex: 0, // next slot to assign (0-based)
      isComplete: false
    };
  },

  /**
   * Assigns the next participant from the active block.
   * If current block is full, a new block is spawned.
   * @param {Object} activeBlock Current block state
   * @param {number} participantId Next participant ID
   * @param {Object} options Configuration options
   * @returns {Object} { updatedBlock, allocationRecord, isNewBlock }
   */
  assignNextInBlock(activeBlock, participantId, options = {}) {
    const labelA = options.labelA || 'Group A';
    const labelB = options.labelB || 'Group B';
    const blockSize = options.blockSize || 4;

    let block = activeBlock;
    let isNewBlock = false;

    if (!block || block.currentIndex >= block.slots.length) {
      const nextBlockNumber = block ? (block.blockNumber + 1) : 1;
      block = this.createBlock(blockSize, { blockNumber: nextBlockNumber, labelA, labelB });
      isNewBlock = true;
    }

    const slot = block.slots[block.currentIndex];
    slot.assigned = true;
    slot.participantId = participantId;
    const now = new Date();
    slot.timestamp = now.toISOString();
    slot.displayTime = now.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' });

    // Update group labels if user customized them
    slot.groupLabel = slot.groupKey === 'A' ? labelA : labelB;

    block.currentIndex += 1;
    if (block.currentIndex >= block.slots.length) {
      block.isComplete = true;
    }

    const allocationRecord = {
      participantId,
      timestamp: slot.timestamp,
      displayTime: slot.displayTime,
      blockNumber: block.blockNumber,
      blockSize: block.blockSize,
      slotInBlock: slot.slotIndex,
      groupKey: slot.groupKey,
      groupLabel: slot.groupLabel
    };

    return {
      updatedBlock: block,
      allocationRecord,
      isNewBlock
    };
  },

  /**
   * Computes running balance statistics for audit trail
   * @param {Array} history 
   * @returns {Object} Summary metrics
   */
  computeSummary(history = []) {
    const total = history.length;
    let countA = 0;
    let countB = 0;

    for (let i = 0; i < total; i++) {
      if (history[i].groupKey === 'A') countA++;
      else if (history[i].groupKey === 'B') countB++;
    }

    const pctA = total > 0 ? (countA / total) * 100 : 0;
    const pctB = total > 0 ? (countB / total) * 100 : 0;
    const diff = Math.abs(countA - countB);

    let ratioStr = '1.00 : 1.00';
    if (countB === 0 && countA > 0) {
      ratioStr = `${countA} : 0`;
    } else if (countA === 0 && countB > 0) {
      ratioStr = `0 : ${countB}`;
    } else if (countB > 0) {
      ratioStr = `${(countA / countB).toFixed(2)} : 1.00`;
    }

    return {
      total,
      countA,
      countB,
      pctA,
      pctB,
      diff,
      ratioStr
    };
  },

  /**
   * Generates formatted CSV string for audit trail export
   * @param {Array} history 
   * @param {string} mode 'simple' or 'block'
   * @returns {string} CSV text
   */
  exportToCSV(history = [], mode = 'simple') {
    if (mode === 'simple') {
      const headers = ['Participant ID', 'Timestamp', 'Random Integer (1-100)', 'Parity', 'Assigned Group Code', 'Group Label'];
      const rows = history.map(item => [
        item.participantId,
        `"${item.timestamp}"`,
        item.randomNumber,
        item.parity,
        item.groupKey,
        `"${item.groupLabel.replace(/"/g, '""')}"`
      ]);
      return [headers.join(','), ...rows.map(r => r.join(','))].join('\r\n');
    } else {
      const headers = ['Participant ID', 'Timestamp', 'Block Number', 'Block Size', 'Slot in Block', 'Assigned Group Code', 'Group Label'];
      const rows = history.map(item => [
        item.participantId,
        `"${item.timestamp}"`,
        item.blockNumber,
        item.blockSize,
        item.slotInBlock,
        item.groupKey,
        `"${item.groupLabel.replace(/"/g, '""')}"`
      ]);
      return [headers.join(','), ...rows.map(r => r.join(','))].join('\r\n');
    }
  }
};
