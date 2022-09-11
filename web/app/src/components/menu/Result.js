import {
  MDBIcon as Icon,
  MDBTable as Table,
  MDBTableBody as TableBody,
  MDBTableHead as TableHead,
  MDBTypography as Typography,
} from "mdbreact";
import { Fragment, useEffect, useState } from "react";
import { Link, useNavigate } from "react-router-dom";
import { usePrimerDriverContext } from "../../contexts/PrimerDriverContext";

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
  }, [PDState.loadedResultsFromModule]);

  useEffect(() => {
    return () => {
      PDDispatch({
        type: "updateLoadedResults",
        payload: false,
      });
    };
  }, []);

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
        <Link to="/" className="btn btn-blue-grey mb-4 mr-3" id="back">
          <Icon fas icon="arrow-left" className="mr-3" />
          main menu
        </Link>
        <Typography tag="h2" className="d-md-inline ml-md-3">
          {(() => {
            if (mode === "CHAR") return "Primer characterization";
            else if (mode === "DNA") return "DNA-based primer design";
            else if (mode === "PRO") return "Protein-based primer design";
            else return null;
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
              {Object.keys(results).map((key, i) => (
                <tr key={i}>
                  <th scope="row">{key}</th>
                  <td>{results[key]["1"]}</td>
                </tr>
              ))}
            </TableBody>
          </Table>
        </div>
      );
    else {
      if (Object.keys(results).length >= 1) {
        return (
          <div>
            <Header />
            <Typography tag="h2" variant="h2-responsive" className="mx-md-2 my-4">
              {Object.keys(results).length} result{Object.keys(results).length !== 1 && "s"}
            </Typography>
            {Object.keys(results).map(
              (item, i) =>
                i < viewLimit && (
                  <Fragment key={i}>
                    <Table responsive className="text-nowrap" size="sm" bordered hover>
                      <TableHead>
                        <tr>
                          <th scope="col" />
                          <th scope="col">{`Primer ${i + 1}`}</th>
                        </tr>
                      </TableHead>
                      <TableBody>
                        {Object.keys(results[(i + 1).toString()]).map((key, j) => (
                          <tr key={j}>
                            <th scope="row">{key}</th>
                            <td>{results[(i + 1).toString()][key]}</td>
                          </tr>
                        ))}
                      </TableBody>
                    </Table>
                    {i + 1 === viewLimit && i + 1 < Object.keys(results).length && (
                      <div className="text-center mt-4 mb-5">
                        <button className="btn btn-info" onClick={handleViewMore}>
                          View more...
                        </button>
                      </div>
                    )}
                  </Fragment>
                ),
            )}
          </div>
        );
      } else {
        return (
          <div>
            <Header />
            <Typography tag="h2" variant="h2-responsive" className="mx-md-2 my-4 text-center">
              {results.data}
            </Typography>
          </div>
        );
      }
    }
  } else {
    return (
      <div className="text-center">
        <Header />
        <p>
          Oops! Something went wrong on the server. Please try again later, or{" "}
          <a href="https://github.com/kvdomingo/primerdriver/issues" target="_blank" rel="noopener noreferrer">
            report an issue
          </a>
          .
        </p>
      </div>
    );
  }
}

export default Result;
