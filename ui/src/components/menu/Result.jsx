import {
  MDBIcon as Icon,
  MDBTable as Table,
  MDBTableBody as TableBody,
  MDBTableHead as TableHead,
  MDBTypography as Typography,
} from "mdbreact";
import { Fragment, useEffect, useState } from "react";
import { Link, useNavigate } from "react-router-dom";
import { usePrimerDriverContext } from "../../contexts/PrimerDriverContext.jsx";

const styles = {
  header: {
    backgroundColor: "white",
    position: "sticky",
    top: 0,
    padding: "1.5em 0 1em 0",
  },
};

function Result() {
  const navigate = useNavigate();
  const [mode, setMode] = useState("");
  const [results, setResults] = useState({});
  const [viewLimit, setViewLimit] = useState(10);
  const { PDState, PDDispatch } = usePrimerDriverContext();

  useEffect(() => {
    if (!PDState.loadedResultsFromModule) navigate("/");
  }, [PDState.loadedResultsFromModule, navigate]);

  useEffect(() => {
    return () => {
      PDDispatch({
        type: "updateLoadedResults",
        payload: false,
      });
    };
  }, [PDDispatch]);

  useEffect(() => {
    if (PDState.results.loaded) {
      setResults(PDState.results.data);
      setMode(PDState.mode);
    }
  }, [PDState.results, PDState.mode]);

  function handleViewMore() {
    setViewLimit(viewLimit + 10);
  }

  function Header() {
    return (
      <div style={styles.header}>
        <Link to="/" className="btn btn-blue-grey mr-3 mb-4" id="back">
          <Icon fas icon="arrow-left" className="mr-3" />
          main menu
        </Link>
        <Typography tag="h2" className="d-md-inline ml-md-3">
          {(() => {
            if (mode === "CHAR") return "Primer characterization";
            if (mode === "DNA") return "DNA-based primer design";
            if (mode === "PRO") return "Protein-based primer design";
            return null;
          })()}
        </Typography>
      </div>
    );
  }

  if (typeof results === "object") {
    if (mode === "CHAR")
      return (
        <div>
          <Header />
          <Table responsive bordered size="sm" hover className="text-nowrap">
            <TableBody>
              {Object.keys(results).map(key => (
                <tr key={key}>
                  <th scope="row">{key}</th>
                  <td>{results[key]["1"]}</td>
                </tr>
              ))}
            </TableBody>
          </Table>
        </div>
      );

    if (Object.keys(results).length >= 1) {
      return (
        <div>
          <Header />
          <Typography tag="h2" variant="h2-responsive" className="mx-md-2 my-4">
            {Object.keys(results).length} result
            {Object.keys(results).length !== 1 && "s"}
          </Typography>
          {Object.keys(results).map(
            (item, i) =>
              i < viewLimit && (
                <Fragment key={item}>
                  <Table responsive className="text-nowrap" size="sm" bordered hover>
                    <TableHead>
                      <tr>
                        <th scope="col" />
                        <th scope="col">{`Primer ${i + 1}`}</th>
                      </tr>
                    </TableHead>
                    <TableBody>
                      {Object.keys(results[(i + 1).toString()]).map(key => (
                        <tr key={key}>
                          <th scope="row">{key}</th>
                          <td>{results[(i + 1).toString()][key]}</td>
                        </tr>
                      ))}
                    </TableBody>
                  </Table>
                  {i + 1 === viewLimit && i + 1 < Object.keys(results).length && (
                    <div className="mt-4 mb-5 text-center">
                      <button
                        className="btn btn-info"
                        onClick={handleViewMore}
                        type="button"
                      >
                        View more...
                      </button>
                    </div>
                  )}
                </Fragment>
              ),
          )}
        </div>
      );
    }

    return (
      <div>
        <Header />
        <Typography
          tag="h2"
          variant="h2-responsive"
          className="mx-md-2 my-4 text-center"
        >
          {results.data}
        </Typography>
      </div>
    );
  }

  return (
    <div className="text-center">
      <Header />
      <p>
        Oops! Something went wrong on the server. Please try again later, or{" "}
        <a
          href="https://github.com/kvdomingo/primerdriver/issues"
          target="_blank"
          rel="noopener noreferrer"
        >
          report an issue
        </a>
        .
      </p>
    </div>
  );
}

export default Result;
