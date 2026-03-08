/**
 * JavaScript port of lemkelcp.R (itself a port of Andy Lamperski's Python impl)
 * Solves LCP: find w, z >= 0 s.t. w - M*z = q, w'*z = 0
 *
 * The R implementation uses column-swapping in the tableau.
 * Columns layout: [w_1..w_n | z_1..z_n | y (driver) | q ]
 * Tind tracks which variable each column represents: [type; index]
 *   type: W=0, Z=1, Y=2, Q=3
 */

const W_TYPE = 0, Z_TYPE = 1, Y_TYPE = 2, Q_TYPE = 3;

function lemkelcp(M, q, maxIter = 200) {
    const n = q.length;

    // Build tableau: n x (2n+2)
    // [I | -M | -1 | q]
    const Tmat = Array.from({ length: n }, (_, i) => {
        const row = new Array(2 * n + 2).fill(0);
        row[i] = 1.0;                          // identity (w cols)
        for (let j = 0; j < n; j++) row[n + j] = -M[i][j]; // -M (z cols)
        row[2 * n] = -1.0;                     // driver y col
        row[2 * n + 1] = q[i];                 // rhs q col
        return row;
    });

    // Tind: 2 x (2n+2), each column is [type, index]
    // First n cols: basic w_i
    // Next n cols:  non-basic z_i
    // col 2n:       driver y
    // col 2n+1:     q (rhs, never swapped)
    const Tind = Array.from({ length: 2 * n + 2 }, (_, c) => {
        if (c < n) return [W_TYPE, c];
        if (c < 2 * n) return [Z_TYPE, c - n];
        if (c === 2 * n) return [Y_TYPE, 0];
        return [Q_TYPE, 0];
    });

    // wPos[i] = current column of w_i
    const wPos = Array.from({ length: n }, (_, i) => i);
    // zPos[i] = current column of z_i
    const zPos = Array.from({ length: n }, (_, i) => n + i);

    const NCOLS = 2 * n + 2;
    const Y_COL_INIT = 2 * n;   // initial column of driver

    // ---------- helpers ----------
    function swapColumns(ci, cj) {
        // swap Tind columns
        const ti = Tind[ci].slice(), tj = Tind[cj].slice();
        Tind[ci] = tj; Tind[cj] = ti;
        // update wPos / zPos
        const updatePos = (type, idx, newCol) => {
            if (type === W_TYPE) wPos[idx] = newCol;
            else if (type === Z_TYPE) zPos[idx] = newCol;
        };
        updatePos(ti[0], ti[1], cj);
        updatePos(tj[0], tj[1], ci);
        // swap Tmat columns
        for (let r = 0; r < n; r++) {
            const tmp = Tmat[r][ci];
            Tmat[r][ci] = Tmat[r][cj];
            Tmat[r][cj] = tmp;
        }
    }

    // clearDriverColumn: pivot row `ind` on column NCOLS-2 (current "entering" col)
    function clearDriverColumn(row) {
        const driverCol = NCOLS - 2;  // column just before q
        const a = Tmat[row][driverCol];
        if (Math.abs(a) < 1e-14) return;
        for (let c = 0; c < NCOLS; c++) Tmat[row][c] /= a;
        for (let r = 0; r < n; r++) {
            if (r !== row) {
                const b = Tmat[r][driverCol];
                for (let c = 0; c < NCOLS; c++) Tmat[r][c] -= b * Tmat[row][c];
            }
        }
    }

    // partnerPos: given a row index into the basis (which variable just left),
    // return the COLUMN of its partner (complementary variable)
    function partnerPos(basisRow) {
        // after pivot, that variable has moved to position NCOLS-2 (was swapped there)
        // But actually we look at Tind of basisRow col = the current "driver" slot... no.
        // In the R code, pivot(pos) looks at ppos = partnerPos(pos) where pos is the ROW
        // of the min-ratio test winner. 
        // But Tind is not indexed by row – it's indexed by column!
        // Wait: in R's impl, the Tind has n+2 columns (basic w, non-basic z, y, q).
        // Each column of Tind says which variable THAT COLUMN represents.
        // After swapColumns, the "driver" is at NCOLS-2.
        // The variable that JUST LEFT is the one in the driver column position?
        // Let me re-read: clearDriverColumn pivots on column NCOLS-2, and is called
        // with the row index `ind`. After that, pivot(ind) does:
        //   ppos = partnerPos(ind)  -> but partnerPos takes a COLUMN position in Tind.
        // In R: partnerPos(pos) looks at Tind[,pos], which is the column `pos` of Tind.
        // The `pos` in init/step is the ROW index of the min-ratio test winner.
        // But Tind columns ≠ rows... 
        // WAIT: In R's Tind, the first n columns correspond to the BASIC variables
        // (currently in the basis). So Tind[,i] = [W, i] means row i of the tableau
        // has w_i as the basic variable. This means Tind columns = tableau rows for basic vars!
        // So partnerPos(pos) with pos = a ROW index looks at what basic var is in that row.
        const [type, idx] = Tind[basisRow]; // column = row (for basic vars)
        if (type === W_TYPE) return zPos[idx];
        if (type === Z_TYPE) return wPos[idx];
        return null; // Y type -> driver leaving
    }

    // pivot: update which variable is "entering" by swapping cols with partner
    function doPivot(row) {
        const ppos = partnerPos(row);
        if (ppos !== null && ppos !== undefined) {
            swapColumns(row, ppos);
            swapColumns(row, NCOLS - 2);
            return true;
        } else {
            // driver (Y) leaves – solution found next check
            swapColumns(row, NCOLS - 2);
            return false;
        }
    }

    // ---------- init ----------
    const qCol = Array.from({ length: n }, (_, i) => Tmat[i][NCOLS - 1]);
    const minQ = Math.min(...qCol);
    if (minQ >= 0) {
        // q >= 0: trivial solution z=0
        return { z: new Array(n).fill(0), code: 0, str: 'Trivial solution' };
    }

    // Find row with most negative q
    let initRow = 0;
    for (let i = 1; i < n; i++) if (qCol[i] < qCol[initRow]) initRow = i;

    clearDriverColumn(initRow);
    doPivot(initRow);

    // ---------- iterate ----------
    for (let k = 1; k <= maxIter; k++) {
        // Check if driver (Y) is no longer basic: Tind at position NCOLS-2 should be Y
        if (Tind[NCOLS - 2][0] === Y_TYPE) {
            // Y is in "non-basic-driver" slot – solution found
            return { z: extractSolution(), code: 0, str: `Found after ${k} iters` };
        }

        // Ratio test on current "driver" column (NCOLS-2)
        const qVals = Array.from({ length: n }, (_, i) => Tmat[i][NCOLS - 1]);
        const aVals = Array.from({ length: n }, (_, i) => Tmat[i][NCOLS - 2]);
        let minRatio = Infinity, pivotRow = -1;
        for (let i = 0; i < n; i++) {
            if (aVals[i] > 1e-12) {
                const ratio = qVals[i] / aVals[i];
                if (ratio < minRatio) { minRatio = ratio; pivotRow = i; }
            }
        }

        if (pivotRow === -1) {
            return { z: null, code: 1, str: 'Secondary ray – no solution' };
        }

        clearDriverColumn(pivotRow);
        doPivot(pivotRow);
    }

    return { z: null, code: 2, str: 'Max iterations exceeded' };

    function extractSolution() {
        const z = new Array(n).fill(0);
        const qVals = Array.from({ length: n }, (_, i) => Tmat[i][NCOLS - 1]);
        for (let i = 0; i < n; i++) {
            if (Tind[i][0] === Z_TYPE) {
                z[Tind[i][1]] = Math.max(0, qVals[i]);
            }
        }
        return z;
    }
}
