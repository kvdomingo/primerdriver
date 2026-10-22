import { Link, useNavigate } from "@tanstack/react-router";
import {
  MDBIcon as Icon,
  MDBTable as Table,
  MDBTableBody as TableBody,
  MDBTableHead as TableHead,
} from "mdbreact";
import { type CSSProperties, Fragment, useEffect, useState } from "react";
import { type ResultData, usePrimerDriverStore } from "../../store";

const styles: Record<string, CSSProperties> = {
  header: {
    backgroundColor: "white",
    position: "sticky",
    top: 0,
    padding: "1.5em 0 1em 0",
  },
};

function isMessage(data: ResultData): data is { data: string } {
  return typeof data === "object" && typeof data.data === "string";
}

function Result() {
  const navigate = useNavigate();
  const [viewLimit, setViewLimit] = useState(10);
  const { loadedResultsFromModule, setLoadedResultsFromModule, ...store } =
    usePrimerDriverStore();
  const results = store.results.loaded ? store.results.data : {};
  const mode = store.results.loaded ? store.mode : "";

  useEffect(() => {
    if (!loadedResultsFromModule) navigate({ to: "/" });
  }, [loadedResultsFromModule, navigate]);

  useEffect(() => {
    return () => setLoadedResultsFromModule(false);
  }, [setLoadedResultsFromModule]);

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
        <h2 className="d-md-inline ml-md-3">
          {(() => {
            if (mode === "CHAR") return "Primer characterization";
            if (mode === "DNA") return "DNA-based primer design";
            if (mode === "PRO") return "Protein-based primer design";
            return null;
          })()}
        </h2>
      </div>
    );
  }

  if (typeof results === "object") {
    if (isMessage(results) || Object.keys(results).length === 0) {
      return (
        <div>
          <Header />
          <h2 className="h2-responsive mx-md-2 my-4 text-center">
            {isMessage(results) ? results.data : null}
          </h2>
        </div>
      );
    }

    if (mode === "CHAR")
      return (
        <div>
          <Header />
          <Table responsive bordered small hover className="text-nowrap">
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

    return (
      <div>
        <Header />
        <h2 className="h2-responsive mx-md-2 my-4">
          {Object.keys(results).length} result
          {Object.keys(results).length !== 1 && "s"}
        </h2>
        {Object.keys(results).map(
          (item, i) =>
            i < viewLimit && (
              <Fragment key={item}>
                <Table responsive className="text-nowrap" small bordered hover>
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
