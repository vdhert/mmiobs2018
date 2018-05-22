pattern = [
            [1.26462, -0.673997, -3.024425],
            [0, 0, -2.5],
            [0, 0, 0],
            [-1.23670, -0.656232, -3.010602]]
        # positions of atoms/points C, C alfa, C beta extended by 1 A and N, respectively
        try:
            coords = [self.atoms[i] for i in ("C", "CA", "CB", "N")]
        except KeyError:
            coords = (self.atoms['C'], self.atoms['CA'], (self.atoms['CA'] - self.atoms['C']) + (self.atoms['CA'] - self.atoms['N']), self.atoms['N'])
        bb_coords = [coord_obj.get_coord() for coord_obj in coords]
        #

        # ===============
        # foregin code starts here
        # ===============
        # all subsequent comments untile notice were made by author
        # assertions replaced by ValueErrors
        # pylint:disable=invalid-name, no-member

        # check for consistency
        if len(bb_coords) != len(pattern):
            raise ValueError('Wrong lenght of backbone: mer %s' % str(self))
        L = len(bb_coords)

        # must alway center the two proteins to avoid
        # affine transformations.  Center the two proteins
        # to their selections.
        COM1 = numpy.sum(bb_coords, axis=0) / float(L)
        COM2 = numpy.sum(pattern, axis=0) / float(L)
        bb_coords = bb_coords - COM1
        pattern = pattern - COM2

        # This beautiful step provides the answer. V and Wt are the orthonormal
        # bases that when multiplied by each other give us the rotation matrix, U.
        # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
        V, S, Wt = numpy.linalg.svd(numpy.dot(numpy.transpose(pattern), bb_coords))

        # we alredy have our solution, in the aaults from SVD.
        # we just need to check for reflections and then produce
        # the rotation.  V and Wt are orthonormal, so their det's
        # are +/-1.
        reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
        if reflect == -1.0:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # U is simply V*Wt
        U = numpy.dot(V, Wt)

        # rotate and translate the molecule
        pattern = numpy.dot((pattern), U) + COM1
        pattern = pattern.tolist()
