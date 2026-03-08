/**
 * Lemke's Algorithm for Linear Complementarity Problem (LCP)
 * Find w, z >= 0 s.t. w - Mz = q and w'z = 0
 */
class LemkeTableau {
    constructor(M, q, maxIter = 100) {
        const n = q.length;
        this.n = n;
        this.maxIter = maxIter;

        // Tmat: [I | -M | -1... | q]
        // Size: n x (2n + 2)
        this.Tmat = Array.from({ length: n }, (_, i) => {
            const row = new Array(2 * n + 2).fill(0);
            row[i] = 1; // I
            for (let j = 0; j < n; j++) {
                row[n + j] = -M[i][j]; // -M
            }
            row[2 * n] = -1; // -1 (driver z0)
            row[2 * n + 1] = q[i]; // q
            return row;
        });

        this.wPos = Array.from({ length: n }, (_, i) => i);
        this.zPos = Array.from({ length: n }, (_, i) => n + i);
        this.W = 0;
        this.Z = 1;
        this.Y = 2; // driver z0
        this.Q = 3;

        // Tind tracks which variable is in each column
        // 0: W, 1: Z, 2: Y, 3: Q
        this.Tind = [];
        // Basic variables indices (0 to n-1)
        this.basis = Array.from({ length: n }, (_, i) => ({ type: this.W, index: i }));
        // Non-basic variables indicators in columns
        this.cols = [];
        for (let i = 0; i < n; i++) this.cols.push({ type: this.W, index: i }); // Basic cols
        for (let i = 0; i < n; i++) this.cols.push({ type: this.Z, index: i }); // Non-basic z
        this.cols.push({ type: this.Y, index: 0 }); // Driver
        this.cols.push({ type: this.Q, index: 0 }); // Q
    }

    // This is a more direct translation of the R6 class logic
    lemkeAlgorithm() {
        if (this.init()) {
            for (let k = 1; k <= this.maxIter; k++) {
                if (!this.step()) {
                    return { z: null, code: 1, str: 'Secondary ray found' };
                }
                if (this.isDriverBasic()) {
                    // This is reached when the "entering" variable was the driver?
                    // Wait, the logic in R6 says col ncol-1 is Y.
                }
                // Check if Y is no longer basic?
                const driverRow = this.basis.findIndex(v => v.type === this.Y);
                if (driverRow === -1) {
                    return { z: this.extractSolution(), code: 0, str: `Solution found after ${k} iterations` };
                }
            }
            return { z: null, code: 2, str: 'Max iterations exceeded' };
        } else {
            return { z: new Array(this.n).fill(0), code: 0, str: 'Solution found after 0 iterations' };
        }
    }

    init() {
        let minQ = 0;
        let minIdx = -1;
        for (let i = 0; i < this.n; i++) {
            if (this.Tmat[i][2 * this.n + 1] < minQ) {
                minQ = this.Tmat[i][2 * this.n + 1];
                minIdx = i;
            }
        }

        if (minIdx !== -1) {
            this.pivot(minIdx, 2 * this.n); // Pivot on driver column (index 2n)
            return true;
        }
        return false;
    }

    step() {
        // Entering variable is the partner of the variable that just left the basis
        const enteringCol = this.nextEnteringCol;

        // Ratio test
        let minRatio = Infinity;
        let leavingRow = -1;
        for (let i = 0; i < this.n; i++) {
            const val = this.Tmat[i][enteringCol];
            if (val > 1e-12) {
                const ratio = this.Tmat[i][2 * this.n + 1] / val;
                if (ratio < minRatio) {
                    minRatio = ratio;
                    leavingRow = i;
                }
            }
        }

        if (leavingRow !== -1) {
            this.pivot(leavingRow, enteringCol);
            return true;
        }
        return false;
    }

    pivot(row, col) {
        const pivotVal = this.Tmat[row][col];
        // Normalize pivot row
        for (let j = 0; j < 2 * this.n + 2; j++) {
            this.Tmat[row][j] /= pivotVal;
        }
        // Eliminate other rows
        for (let i = 0; i < this.n; i++) {
            if (i !== row) {
                const val = this.Tmat[i][col];
                for (let j = 0; j < 2 * this.n + 2; j++) {
                    this.Tmat[i][j] -= val * this.Tmat[row][j];
                }
            }
        }

        // Update basis
        const leavingVar = this.basis[row];
        const enteringVar = this.cols[col];
        this.basis[row] = enteringVar;

        // Determine next entering column
        // If w_i left, z_i enters. If z_i left, w_i enters.
        if (leavingVar.type === this.W) {
            this.nextEnteringCol = this.n + leavingVar.index;
        } else if (leavingVar.type === this.Z) {
            this.nextEnteringCol = leavingVar.index;
        } else {
            // Driver left! Solution found.
            this.nextEnteringCol = -1;
        }
    }

    extractSolution() {
        const z = new Array(this.n).fill(0);
        for (let i = 0; i < this.n; i++) {
            if (this.basis[i].type === this.Z) {
                z[this.basis[i].index] = Math.max(0, this.Tmat[i][2 * this.n + 1]);
            }
        }
        return z;
    }
}

function lemkelcp(M, q, maxIter = 100) {
    const tableau = new LemkeTableau(M, q, maxIter);
    return tableau.lemkeAlgorithm();
}
